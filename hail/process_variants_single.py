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
    # 增加 encoding='ISO-8859-1' 解决 UnicodeDecodeError
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = [c.strip() for c in df.columns]
    
    idx = pd.MultiIndex.from_frame(df[key_cols])
    df["__vals__"] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df["__vals__"].values, index=idx).to_dict()


def load_simple_val_db(path, key_cols, val_col, delimiter: str = "\t") -> dict:
    """Load a tabular file as a dict returning a single value."""
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = [c.strip() for c in df.columns]
    
    missing = [c for c in key_cols + [val_col] if c not in df.columns]
    if missing:
        print(f"ERROR: Missing columns {missing} in {path}")
        print(f"Available columns: {list(df.columns)}")
        raise KeyError(f"Columns {missing} not found in {path}")

    idx = pd.MultiIndex.from_frame(df[key_cols])
    return pd.Series(df[val_col].values, index=idx).to_dict()


def load_pos_only_db(path, pos_col, val_col, delimiter: str = "\t") -> dict:
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = df.columns.str.strip()
    return pd.Series(df[val_col].values, index=df[pos_col]).to_dict()


# ==============================================================================
# Core
# ==============================================================================

def process_and_filter_variants(args) -> None:
    print("\n--- Starting single-sample variant processing ---")

    min_vaf = float(args.min_vaf) if args.min_vaf else 0.01
    min_dp = float(args.min_dp) if args.min_dp else 100
    do_post_filter = str(args.post_filtering).lower() == "true"

    # ---- metadata (haplogroup, contamination) ----
    hap_dict = pd.read_csv(
        args.fullhaplogroups, sep="\t", header=0, dtype=str, index_col=0
    ).squeeze("columns").to_dict()

    contam_dict = pd.read_csv(
        args.contamination, sep="\t", header=0, dtype=str, index_col=0
    ).squeeze("columns").to_dict()

    # Optional: disease meta
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

    hgv = pd.read_csv(args.haplogroup_varcache, sep="\t", dtype=str)
    haplovar_db = pd.Series(hgv.Assoc_haplogroups.values, index=hgv.Variant) \
                    .str.lower().to_dict()

    apogee_db = load_complex_db(
        args.mitimpactcache, ["Ref", "Start", "Alt"], ["APOGEE1", "APOGEE2"]
    )

    # 修改：根据您的文件头，使用大写 'Ref', 'Start', 'Alt'
    napogee_db = load_simple_val_db(args.napogeecache, ["ref", "start", "alt"], "nAPOGEE_score")
    tapogee_db = load_simple_val_db(args.tapogeecache, ["Ref", "Pos", "Alt"], "t-APOGEE score")

    # 修改：MLC 数据库加载
    # SNV 使用完全匹配
    mlc_snv_db = load_simple_val_db(args.mlc_snv_cache, ["Reference", "Position", "Alternate"], "MLC_score")
    # Indel 仅使用位置匹配 (假设 Indel 文件中位置列名为 'Position')
    mlc_indel_db = load_pos_only_db(args.mlc_indel_cache, "Position", "MLC_pos_score")

    # NEW: Regional Constraint
    rc_df = pd.read_csv(args.regional_constraint_cache, sep="\t", dtype=str, encoding='ISO-8859-1')
    rc_df.columns = [c.strip() for c in rc_df.columns]
    rc_db = pd.Series(
        list(zip(rc_df["in_rc"], rc_df["min_distance_to_rc"])),
        index=pd.MultiIndex.from_frame(rc_df[["REF", "POS", "ALT"]])
    ).to_dict()

    hmtvar_db = load_complex_db(
        args.hmtvarcache, ["REF", "POS", "ALT"], ["HmtVar"]
    )

    # MitoTIP
    mitotip_map = {"Q1": "likely pathogenic", "Q2": "possibly pathogenic", "Q3": "possibly benign", "Q4": "likely benign"}
    mitotip_df = pd.read_csv(args.mitotipcache, sep="\t", dtype=str)
    mitotip_df["prediction"] = mitotip_df["Quartile"].map(mitotip_map)
    mitotip_db = pd.Series(
        mitotip_df["prediction"].values,
        index=tuple(zip(mitotip_df["rCRS"], mitotip_df["Position"], mitotip_df["Alt"]))
    ).to_dict()

    # HelixMTdb
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

    # ClinVar
    clin = pd.read_csv(args.clinvarcache, sep="\t", dtype=str, encoding='ISO-8859-1')
    clin.columns = [c.strip() for c in clin.columns]
    
    gene_col = "Gene(s)" if "Gene(s)" in clin.columns else "Gene"
    
    clin["pos"] = clin["GRCh38Location"]
    spdi = clin["Canonical SPDI"].str.split(":", expand=True)
    clin["ref"], clin["alt"] = spdi[2], spdi[3]
    
    clin[gene_col] = clin[gene_col].fillna("").str.strip()

    clin = clin[
        (clin["ref"].str.len() == 1) &
        (clin["alt"].str.len() == 1) &
        (clin["ref"] != clin["alt"]) &
        (clin["Germline classification"] != "Conflicting interpretations of pathogenicity")
    ].copy()

    clinvar_db = pd.Series(
        clin["Germline classification"].values,
        index=pd.MultiIndex.from_frame(clin[["ref", "pos", "alt", gene_col]])
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

    sample_col = vcf.columns[9]
    vcf = vcf.rename(columns={sample_col: "SAMPLE_DATA"})

    if str(contam_dict.get(sample_col, "")).upper() == "YES":
        print(f"[!] Sample {sample_col} flagged as contaminated; no output.")
        return

    vcf = vcf[~vcf["FILTER"].str.contains(
        "base_qual|strand_bias|weak_evidence|blacklisted_site|contamination|position",
        na=False
    )]
    if vcf.empty:
        print("[!] No variants after FILTER cleanup.")
        return

    fmt = vcf["SAMPLE_DATA"].str.split(":", n=4, expand=True)
    fmt.columns = ["GT", "AD", "AF", "DP", "Other"][:fmt.shape[1]]
    for c in fmt.columns:
        vcf[c] = fmt[c]
    vcf["Heteroplasmy"] = pd.to_numeric(vcf["AF"], errors="coerce")
    vcf["DP"] = pd.to_numeric(vcf["DP"], errors="coerce")
    vcf = vcf[
        (vcf["Heteroplasmy"].fillna(0) >= min_vaf) & 
        (vcf["DP"].fillna(0) >= min_dp)
    ].copy()
    
    if vcf.empty:
        print(f"[!] No variants passed the hard filters (VAF >= {min_vaf}, DP >= {min_dp}).")
        return
    
    artifact_sites = ["301", "302", "310", "316", "3107", "16182"]
    vcf["is_artifact_prone"] = vcf["POS"].astype(str).isin(artifact_sites)

    print("[*] Calculating indel_stack metrics...")
    
    vcf["is_indel"] = vcf["REF"].str.len() != vcf["ALT"].str.len()
    
    unique_indels = vcf[vcf["is_indel"]][["POS", "ALT"]].drop_duplicates()
    pos_indel_counts = unique_indels.groupby("POS")["ALT"].nunique()
    
    stack_positions = pos_indel_counts[pos_indel_counts >= 2].index.tolist()
    
    vcf["indel_stack"] = vcf["is_indel"] & vcf["POS"].isin(stack_positions)
    
    def update_filter_str(row):
        current_filter = str(row["FILTER"])
        if row["indel_stack"]:
            if current_filter == "PASS" or current_filter == ".":
                return "indel_stack"
            else:
                return current_filter + ";indel_stack"
        return current_filter

    vcf["FILTER"] = vcf.apply(update_filter_str, axis=1)

    def extract_csq_list(info_str):
        for segment in info_str.split(";"):
            if segment.startswith("CSQ="):
                return segment[4:].split(",")
        return [""]

    vcf["CSQ_LIST"] = vcf["INFO"].apply(extract_csq_list)
    
    vcf = vcf.explode("CSQ_LIST").reset_index(drop=True)

    info_cols = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
        "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
        "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "HGVS_OFFSET",
    ]
    
    vep_split = vcf["CSQ_LIST"].str.split("|", expand=True)
    vep_start = 0 
    take = min(len(info_cols), max(0, vep_split.shape[1] - vep_start))
    
    if take > 0:
        vcf[info_cols[:take]] = vep_split.iloc[:, vep_start: vep_start + take]

    vcf.drop(columns=["CSQ_LIST"], inplace=True)

    vcf["SAMPLE_ID"] = sample_col
    vcf["Haplogroup"] = vcf["SAMPLE_ID"].map(hap_dict)

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

    # New Annotation Mapping
    vcf["nAPOGEE"] = [napogee_db.get(k, "") for k in keys]
    vcf["tAPOGEE"] = [tapogee_db.get(k, "") for k in keys]

    def get_mlc_score(row):
        # 如果是 SNV，使用 REF+POS+ALT 匹配
        if row["VARIANT_CLASS"] == "SNV":
            key = (row["REF"], row["POS"], row["ALT"])
            return mlc_snv_db.get(key, "")
        else:
            # 如果是 Indel，仅根据 POS 匹配
            return mlc_indel_db.get(row["POS"], "")

    vcf["MLC_score"] = vcf.apply(get_mlc_score, axis=1)

    # Regional Constraint logic
    def get_rc_status(row):
        cons = str(row.get("Consequence", "")).lower()
        if "missense_variant" not in cons and "rrna_variant" not in cons:
            return "NA"
        
        key = (row["REF"], row["POS"], row["ALT"])
        rc_info = rc_db.get(key)
        if not rc_info:
            return "no"
            
        in_rc = str(rc_info[0]).lower()
        try:
            dist = float(rc_info[1])
        except (ValueError, TypeError):
            dist = 999.0
            
        if in_rc == "yes":
            return "yes"
        elif in_rc == "no" and dist <= 6:
            return "proximal"
        else:
            return "no"

    vcf["in_regional_constraint"] = vcf.apply(get_rc_status, axis=1)

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

    for c in ["gnomad_af_hom", "helix_af_hom", "mitomap_af"]:
        vcf[c] = pd.to_numeric(vcf[c], errors="coerce")

    # Output columns
    final_cols = [
        "SAMPLE_ID", "CHROM", "POS", "REF", "ALT", "FILTER",
        "GT", "AD", "Heteroplasmy", "DP",
        "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "HGVSp",
        "Codons", "VARIANT_CLASS", "indel_stack","is_artifact_prone",
        "Haplogroup", "Haplogroup_Var_Status",
        "gnomad_max_hl", "gnomad_af_hom", "gnomad_af_het",
        "apogee_class", "mitotip_class", "hmtvar_class",
        "nAPOGEE", "tAPOGEE", "MLC_score", "in_regional_constraint",
        "helix_max_hl", "helix_af_hom", "helix_af_het",
        "mitomap_gbcnt", "mitomap_af",
        "mitomap_status", "mitomap_plasmy", "mitomap_disease",
        "clinvar_interp",
    ]
    for c in final_cols:
        if c not in vcf.columns:
            vcf[c] = pd.NA

    os.makedirs(args.final_output_dir, exist_ok=True)

    prefilter_mask = (vcf["FILTER"].isin(["PASS", ".", ""])) & (~vcf["indel_stack"])
    vcf_pre = vcf[prefilter_mask].copy()

    os.makedirs(args.final_output_dir, exist_ok=True)
    pre_path = os.path.join(args.final_output_dir, "variant_list_prefiltering.txt")
    
    vcf_pre[final_cols].to_csv(pre_path, sep="\t", index=False, na_rep="")
    print(f"[+] Prefiltering table (PASS only) saved to: {pre_path}")

    if do_post_filter:
        filtered = vcf[
            (vcf["gnomad_af_hom"] < 0.01) &
            (vcf["helix_af_hom"]  < 0.01) &
            (vcf["mitomap_af"]    < 0.01) &
            (vcf["Consequence"]   != "synonymous_variant") &
            (vcf["Heteroplasmy"].fillna(0) > 0.05) &
            (~vcf["clinvar_interp"].isin(["Benign", "Likely benign", "Affects"])) &
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
    ap.add_argument("--napogeecache", required=True)
    ap.add_argument("--tapogeecache", required=True)
    ap.add_argument("--mlc_snv_cache", required=True)
    ap.add_argument("--mlc_indel_cache", required=True)
    ap.add_argument("--regional_constraint_cache", required=True)
    ap.add_argument("--mitotipcache", required=True)
    ap.add_argument("--hmtvarcache", required=True)
    ap.add_argument("--min-vaf", help="Minimum VAF threshold", default="0.01")
    ap.add_argument("--min-dp", help="Minimum DP threshold", default="100")
    ap.add_argument("--post_filtering", help="Post filtering flag", default="False")
    args = ap.parse_args()

    process_and_filter_variants(args)