#!/usr/bin/env python3

import os
import glob
import io
import json
import argparse
import pandas as pd


# ==============================================================================
# Utilities
# ==============================================================================

def load_complex_db(path, key_cols, val_cols, delimiter: str = "\t") -> dict:
    """Load a tabular annotation file as a dict keyed by multiple columns."""
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#")
    idx = pd.MultiIndex.from_frame(df[key_cols])
    df["__vals__"] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df["__vals__"].values, index=idx).to_dict()


# ==============================================================================
# Core
# ==============================================================================

def process_and_filter_variants(args) -> None:
    print("\n--- Starting single-sample variant processing ---")

    # ---- metadata (haplogroup, contamination) ----
    hap_dict = pd.read_csv(
        args.fullhaplogroups, sep="\t", header=0, dtype=str, index_col=0
    ).squeeze("columns").to_dict()

    contam_dict = pd.read_csv(
        args.contamination, sep="\t", header=0, dtype=str, index_col=0
    ).squeeze("columns").to_dict()

    # Optional: disease meta (SampleID/Category with flexible casing)
    category_dict = {}
    if args.disease_meta_file and os.path.isfile(args.disease_meta_file):
        meta = pd.read_csv(args.disease_meta_file, sep="\t", dtype=str)
        meta.columns = [c.strip().lower() for c in meta.columns]
        if {"sampleid", "category"}.issubset(meta.columns):
            category_dict = pd.Series(
                meta["category"].values, index=meta["sampleid"]
            ).to_dict()

    # ---- external annotation DBs ----
    gnomad_db = load_complex_db(
        args.gnomadcache, ["ref", "position", "alt"],
        ["max_observed_heteroplasmy", "AF_hom", "AF_het"]
    )

    mitomap_poly_db = load_complex_db(
        args.mitomap_polycache, ["ref", "pos", "alt"], ["gbcnt"]
    )

    mitomap_disease_db = load_complex_db(
        args.mitomap_diseasecache, ["ref", "pos", "alt"],
        ["status", "homoplasmy", "heteroplasmy", "disease"]
    )

    # Haplogroup-specific variants (string key REF+POS+ALT → comma-separated haplos)
    hgv = pd.read_csv(args.haplogroup_varcache, sep="\t", dtype=str)
    haplovar_db = pd.Series(hgv.Assoc_haplogroups.values, index=hgv.Variant) \
                    .str.lower().to_dict()

    apogee_db = load_complex_db(
        args.mitimpactcache, ["Ref", "Start", "Alt"], ["APOGEE1", "APOGEE2"]
    )

    hmtvar_db = load_complex_db(
        args.hmtvarcache, ["REF", "POS", "ALT"], ["HmtVar"]
    )

    # MitoTIP (quartile → label)
    mitotip_map = {
        "Q1": "likely pathogenic",
        "Q2": "possibly pathogenic",
        "Q3": "possibly benign",
        "Q4": "likely benign",
    }
    mitotip_df = pd.read_csv(args.mitotipcache, sep="\t", dtype=str)
    mitotip_df["prediction"] = mitotip_df["Quartile"].map(mitotip_map)
    mitotip_db = pd.Series(
        mitotip_df["prediction"].values,
        index=tuple(zip(mitotip_df["rCRS"], mitotip_df["Position"], mitotip_df["Alt"]))
    ).to_dict()

    # HelixMTdb (normalize fields, derive max heteroplasmy)
    helix = pd.read_csv(args.helixcache, sep="\t", dtype=str)
    helix = helix[helix["alleles"].str.count(",") == 1].copy()
    helix["pos"] = helix["locus"].str.split("chrM:").str[1]
    sp = helix["alleles"].str.split('"').str
    helix["ref"], helix["alt"] = sp[1], sp[3]
    for c in ["AF_hom", "AF_het", "max_ARF"]:
        helix[c] = pd.to_numeric(helix.get(c, 0), errors="coerce")
    helix["max_het"] = helix.apply(
        lambda r: 1.0 if (pd.notna(r["AF_hom"]) and r["AF_hom"] > 0)
        else (r["max_ARF"] if pd.notna(r["max_ARF"]) else 0.0),
        axis=1,
    )
    helix_db = pd.Series(
        list(zip(helix["max_het"], helix["AF_hom"], helix["AF_het"])),
        index=pd.MultiIndex.from_frame(helix[["ref", "pos", "alt"]])
    ).to_dict()

    # ClinVar (bi-allelic SNVs, exclude “Conflicting”)
    clin = pd.read_csv(args.clinvarcache, sep="\t", dtype=str)
    clin["pos"] = clin["GRCh38Location"]
    spdi = clin["Canonical SPDI"].str.split(":", expand=True)
    clin["ref"], clin["alt"] = spdi[2], spdi[3]
    clin = clin[
        (clin["ref"].str.len() == 1) &
        (clin["alt"].str.len() == 1) &
        (clin["ref"] != clin["alt"]) &
        (clin["Germline classification"] != "Conflicting interpretations of pathogenicity")
    ].copy()
    clinvar_db = pd.Series(
        clin["Germline classification"].values,
        index=pd.MultiIndex.from_frame(clin[["ref", "pos", "alt"]])
    ).to_dict()

    # ---- read one VEP VCF ----
    vcf_list = sorted(glob.glob(os.path.join(args.vep_vcf_dir, "*_vep.vcf")))
    if not vcf_list:
        print(f"[!] No *_vep.vcf found in {args.vep_vcf_dir}. Nothing to do.")
        return

    vcf_path = vcf_list[0]
    print(f"[*] Using VEP VCF: {vcf_path}")

    with open(vcf_path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    vcf = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str) \
             .rename(columns={"#CHROM": "CHROM"})
    if vcf.empty:
        print("[!] Empty VCF after header removal.")
        return

    # Single-sample: take the first sample column and rename to SAMPLE_DATA
    sample_col = vcf.columns[9]
    vcf = vcf.rename(columns={sample_col: "SAMPLE_DATA"})

    # Hard drop if contaminated
    if str(contam_dict.get(sample_col, "")).upper() == "YES":
        print(f"[!] Sample {sample_col} flagged as contaminated; no output.")
        return

    # Drop upstream hard-fail sites
    vcf = vcf[~vcf["FILTER"].str.contains(
        "base_qual|strand_bias|weak_evidence|blacklisted_site|contamination|position",
        na=False
    )]
    if vcf.empty:
        print("[!] No variants after FILTER cleanup.")
        return

    # FORMAT parsing: “AF” → Heteroplasmy; DP to numeric
    fmt = vcf["SAMPLE_DATA"].str.split(":", n=4, expand=True)
    fmt.columns = ["GT", "AD", "AF", "DP", "Other"][:fmt.shape[1]]
    for c in fmt.columns:
        vcf[c] = fmt[c]
    vcf["Heteroplasmy"] = pd.to_numeric(vcf["AF"], errors="coerce")
    vcf["DP"] = pd.to_numeric(vcf["DP"], errors="coerce")

    # ---- FIXED VEP INFO parsing ----
    # VCF INFO fields often contain tags like DP=X;AF=Y;CSQ=Allele|Consequence...
    # Direct splitting of the whole INFO column causes index shifting.
    # We must isolate the CSQ (VEP) string first.
    def extract_csq(info_str):
        for segment in info_str.split(";"):
            if segment.startswith("CSQ="):
                return segment[4:] # Remove 'CSQ=' prefix
        return ""

    vcf["CSQ_STRING"] = vcf["INFO"].apply(extract_csq)
    
    info_cols = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
        "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
        "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "HGVS_OFFSET",
    ]
    
    # Split the isolated CSQ string by the pipe character
    vep_split = vcf["CSQ_STRING"].str.split("|", expand=True)
    
    # After stripping 'CSQ=', index 0 is Allele, 1 is Consequence, 2 is IMPACT
    vep_start = 0 
    take = min(len(info_cols), max(0, vep_split.shape[1] - vep_start))
    
    if take > 0:
        vcf[info_cols[:take]] = vep_split.iloc[:, vep_start: vep_start + take]

    # Clean up temporary parsing column
    vcf.drop(columns=["CSQ_STRING"], inplace=True)

    # Metadata columns
    vcf["SAMPLE_ID"] = sample_col
    vcf["Haplogroup"] = vcf["SAMPLE_ID"].map(hap_dict)

    # Disease meta (if present)
    if category_dict:
        vcf["Sample_Category"] = vcf["SAMPLE_ID"].str.split("-").str[0] \
            .str.lower().map(category_dict)
    else:
        vcf["Sample_Category"] = pd.NA

    # Annotation lookups
    keys = list(zip(vcf["REF"], vcf["POS"], vcf["ALT"]))
    vcf["gnomad_max_hl"], vcf["gnomad_af_hom"], vcf["gnomad_af_het"] = zip(
        *[gnomad_db.get(k, ("0", "0", "0")) for k in keys]
    )
    vcf["helix_max_hl"], vcf["helix_af_hom"], vcf["helix_af_het"] = zip(
        *[helix_db.get(k, (0.0, 0.0, 0.0)) for k in keys]
    )
    vcf["clinvar_interp"] = [clinvar_db.get(k, "") for k in keys]
    vcf["mitomap_gbcnt"] = [mitomap_poly_db.get(k, ("0",))[0] for k in keys]
    vcf["mitomap_af"] = pd.to_numeric(vcf["mitomap_gbcnt"], errors="coerce") / 61134

    md = [mitomap_disease_db.get(k, ("", "", "", "")) for k in keys]
    vcf["mitomap_status"] = [d[0] for d in md]
    vcf["mitomap_plasmy"] = [f"{d[1]}/{d[2]}" if d[1] or d[2] else "" for d in md]
    vcf["mitomap_disease"] = [d[3] for d in md]

    # Haplogroup-variant status (presence vs. current sample’s haplogroup)
    def haplo_status(row):
        key = f"{row['REF']}{row['POS']}{row['ALT']}"
        assoc = haplovar_db.get(key)
        if assoc:
            return (
                "haplo_var_match"
                if str(row.get("Haplogroup", "")).lower() in assoc
                else "haplo_var_diff_haplo"
            )
        return "not_haplo_var"

    vcf["Haplogroup_Var_Status"] = vcf.apply(haplo_status, axis=1)

    # APOGEE / MitoTIP / HmtVar
    vcf["apogee_class"] = [
        str(apogee_db.get(k, "")).strip("()").replace("'", "").replace(", ", "/")
        for k in keys
    ]
    vcf["mitotip_class"] = [mitotip_db.get(k, "") for k in keys]

    def get_hmtvar(k):
        val = hmtvar_db.get(k)
        if val and val[0]:
            try:
                return json.loads(val[0]).get("pathogenicity", "")
            except json.JSONDecodeError:
                return ""
        return ""

    vcf["hmtvar_class"] = [get_hmtvar(k) for k in keys]

    # Numeric coercions for filtering/sorting
    for c in ["gnomad_af_hom", "helix_af_hom", "mitomap_af"]:
        vcf[c] = pd.to_numeric(vcf[c], errors="coerce")

    # Single-sample compatibility column
    if "Freq" not in vcf.columns:
        vcf["Freq"] = 0

    # ---- output columns ----
    final_cols = [
        "SAMPLE_ID", "CHROM", "POS", "REF", "ALT", "FILTER",
        "GT", "AD", "Heteroplasmy", "DP",
        "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "HGVSp",
        "Codons", "VARIANT_CLASS",
        "Haplogroup", "Haplogroup_Var_Status",
        "gnomad_max_hl", "gnomad_af_hom", "gnomad_af_het",
        "apogee_class", "mitotip_class", "hmtvar_class",
        "helix_max_hl", "helix_af_hom", "helix_af_het",
        "mitomap_gbcnt", "mitomap_af",
        "mitomap_status", "mitomap_plasmy", "mitomap_disease",
        "clinvar_interp", "Freq",
    ]
    for c in final_cols:
        if c not in vcf.columns:
            vcf[c] = pd.NA

    os.makedirs(args.final_output_dir, exist_ok=True)

    # Prefilter table
    pre_path = os.path.join(args.final_output_dir, "variant_list_prefiltering.txt")
    vcf[final_cols].to_csv(pre_path, sep="\t", index=False, na_rep="")
    print(f"[+] Prefiltering table saved to: {pre_path}")

    # Final filtered table
    filtered = vcf[
        (vcf["gnomad_af_hom"] < 0.01) &
        (vcf["helix_af_hom"]  < 0.01) &
        (vcf["mitomap_af"]    < 0.01) &
        (vcf["Consequence"]   != "synonymous_variant") &
        (pd.to_numeric(vcf["Freq"], errors="coerce").fillna(0) < 100) &
        (vcf["Heteroplasmy"].fillna(0) > 0.05) &
        (~vcf["clinvar_interp"].isin(["Benign", "Likely benign"])) &
        (vcf["Haplogroup_Var_Status"] != "haplo_var_match")
    ].copy()

    filtered = filtered.sort_values(
        by=["gnomad_af_hom", "helix_af_hom", "mitomap_af"],
        ascending=[True, True, True],
        na_position="last"
    )

    out_path = os.path.join(args.final_output_dir, "variant_list.txt")
    filtered[final_cols].to_csv(out_path, sep="\t", index=False, na_rep="")
    print(f"[+] Single-sample variant list saved to: {out_path}")


# ==============================================================================
# CLI
# ==============================================================================

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Single-sample processing of a VEP-annotated VCF.")
    ap.add_argument("--vep-vcf-dir", required=True)
    ap.add_argument("--final-output-dir", required=True)
    ap.add_argument("--fullhaplogroups", required=True)
    ap.add_argument("--contamination", required=True)
    ap.add_argument("--disease-meta-file", required=False, default=None)
    ap.add_argument("--gnomadcache", required=True)
    ap.add_argument("--clinvarcache", required=True)
    ap.add_argument("--mitomap-polycache", required=True)
    ap.add_argument("--mitomap-diseasecache", required=True)
    ap.add_argument("--helixcache", required=True)
    ap.add_argument("--haplogroup-varcache", required=True)
    ap.add_argument("--mitimpactcache", required=True)
    ap.add_argument("--mitotipcache", required=True)
    ap.add_argument("--hmtvarcache", required=True)
    args = ap.parse_args()

    process_and_filter_variants(args)