#!/usr/bin/env python3

import os
import io
import json
import argparse
import pandas as pd
import numpy as np

# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

def load_complex_db(path: str, key_cols: list[str], val_cols: list[str], delimiter: str = "\t") -> dict:
    """Loads external annotation databases into a dictionary for fast lookup."""
    if not os.path.exists(path):
        print(f"WARNING: Database {path} not found. Lookups will return default values.")
        return {}
    
    # 增加 encoding='ISO-8859-1' 解决 UnicodeDecodeError
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    
    # 清洗列名，去除不可见空格
    df.columns = df.columns.str.strip()
    
    missing = [c for c in key_cols + val_cols if c not in df.columns]
    if missing:
        print(f"ERROR: {os.path.basename(path)} missing columns: {missing}")
        print(f"Available columns: {list(df.columns)}")
        raise SystemExit(1)
        
    idx = pd.MultiIndex.from_frame(df[key_cols])
    df["__values__"] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df["__values__"].values, index=idx).to_dict()

def load_pos_only_db(path: str, pos_col: str, val_col: str, delimiter: str = "\t") -> dict:
    """专门为 MLC Indel 设计：仅根据位置加载评分字典"""
    if not os.path.exists(path):
        print(f"WARNING: Database {path} not found.")
        return {}
    
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = df.columns.str.strip()
    
    if pos_col not in df.columns or val_col not in df.columns:
        print(f"ERROR: {os.path.basename(path)} missing columns: {[pos_col, val_col]}")
        raise SystemExit(1)
        
    return pd.Series(df[val_col].values, index=df[pos_col]).to_dict()

def get_hmtvar(k, hmtvar_db):
    """Extracts pathogenicity prediction from HmtVar JSON-like strings."""
    val = hmtvar_db.get(k)
    if val and val[0]:
        try:
            return json.loads(val[0]).get("pathogenicity", "")
        except (json.JSONDecodeError, TypeError):
            return ""
    return ""

# --------------------------------------------------------------------------- #
# Main Processing Logic                                                       #
# --------------------------------------------------------------------------- #

def main() -> None:
    ap = argparse.ArgumentParser(description="Process a VEP-annotated multi-sample VCF and generate summary tables.")
    ap.add_argument("--vep-vcf", required=True, help="Path to one VEP-annotated VCF.")
    ap.add_argument("--final-output-dir", required=True, help="Directory for output tables.")

    # External database paths
    ap.add_argument("--gnomadcache", required=True)
    ap.add_argument("--clinvarcache", required=True)
    ap.add_argument("--mitomap-polycache", required=True)
    ap.add_argument("--mitomap-diseasecache", required=True)
    ap.add_argument("--helixcache", required=True)
    ap.add_argument("--mitimpactcache", required=True)
    ap.add_argument("--mitotipcache", required=True)
    ap.add_argument("--hmtvarcache", required=True)
    
    # NEW arguments
    ap.add_argument("--napogeecache", required=True)
    ap.add_argument("--tapogeecache", required=True)
    ap.add_argument("--mlc_snv_cache", required=True)
    ap.add_argument("--mlc_indel_cache", required=True)
    ap.add_argument("--regional_constraint_cache", required=True)

    # Execution Mode
    ap.add_argument("--pipeline-mode", choices=["population", "disease"], required=True)
    ap.add_argument(
        "--disease-meta-file",
        help="TSV with SampleID and Category columns; Category=Proband rows are kept in disease mode."
    )

    args = ap.parse_args()
    outdir = args.final_output_dir
    os.makedirs(outdir, exist_ok=True)

    # -------------------- Step 1: Handle Disease Mode Meta ----------------- #
    proband_set: set[str] | None = None
    if args.pipeline_mode == "disease":
        if not args.disease_meta_file:
            raise SystemExit("ERROR: --pipeline-mode=disease requires --disease-meta-file.")
        meta = pd.read_csv(args.disease_meta_file, sep="\t", dtype=str)
        if meta.empty:
            raise SystemExit("ERROR: disease_meta_file is empty.")
        meta.columns = [c.strip() for c in meta.columns]
        if not {"SampleID", "Category"}.issubset(meta.columns):
            raise SystemExit("ERROR: disease_meta_file must contain SampleID and Category.")
        proband_set = set(
            meta.loc[meta["Category"].astype(str).str.lower() == "proband", "SampleID"].astype(str)
        )

    # -------------------- Step 2: Load Annotation Databases ---------------- #
    gnomad_db          = load_complex_db(args.gnomadcache,          ["ref", "position", "alt"], ["max_observed_heteroplasmy", "AF_hom", "AF_het"])
    mitomap_poly_db    = load_complex_db(args.mitomap_polycache,    ["ref", "pos", "alt"],      ["gbcnt"])
    mitomap_disease_db = load_complex_db(args.mitomap_diseasecache, ["ref", "pos", "alt"],      ["status", "homoplasmy", "heteroplasmy", "disease"])
    apogee_db          = load_complex_db(args.mitimpactcache,       ["Ref", "Start", "Alt"],    ["APOGEE1", "APOGEE2"])
    hmtvar_db          = load_complex_db(args.hmtvarcache,          ["REF", "POS", "ALT"],      ["HmtVar"])
    
    # NEW: nAPOGEE, tAPOGEE, MLC, and Regional Constraint
    napogee_db = load_complex_db(args.napogeecache, ["ref", "start", "alt"], ["nAPOGEE_score"])
    tapogee_db = load_complex_db(args.tapogeecache, ["Ref", "Pos", "Alt"], ["t-APOGEE score"])
    
    # MLC SNV 使用完全匹配 (假设列名为 Reference, Position, Alternate)
    mlc_snv_db = load_complex_db(args.mlc_snv_cache, ["Reference", "Position", "Alternate"], ["MLC_score"])

    # MLC Indel 仅使用 Position 匹配 (假设 Indel 文件中列名为 Position 和 MLC_pos_score)
    mlc_indel_db = load_pos_only_db(args.mlc_indel_cache, "Position", "MLC_pos_score")
    
    rc_df = pd.read_csv(args.regional_constraint_cache, sep="\t", dtype=str, encoding='ISO-8859-1')
    rc_df.columns = [c.strip() for c in rc_df.columns]
    rc_db = pd.Series(
        list(zip(rc_df["in_rc"], rc_df["min_distance_to_rc"])),
        index=pd.MultiIndex.from_frame(rc_df[["REF", "POS", "ALT"]])
    ).to_dict()

    # Custom Parsing for ClinVar
    clinvar_df = pd.read_csv(args.clinvarcache, sep="\t", dtype=str, encoding='ISO-8859-1')
    clinvar_df.columns = [c.strip() for c in clinvar_df.columns]
    clinvar_df["pos"] = clinvar_df["GRCh38Location"]
    spdi = clinvar_df["Canonical SPDI"].str.split(":", expand=True)
    clinvar_df["ref"], clinvar_df["alt"] = spdi[2], spdi[3]
    clinvar_df = clinvar_df[
        (clinvar_df["ref"].str.len() == 1) & (clinvar_df["alt"].str.len() == 1) &
        (clinvar_df["ref"] != clinvar_df["alt"]) &
        (clinvar_df["Germline classification"] != "Conflicting interpretations of pathogenicity")
    ].copy()
    clinvar_db = pd.Series(
        clinvar_df["Germline classification"].values,
        index=pd.MultiIndex.from_frame(clinvar_df[["ref", "pos", "alt"]])
    ).to_dict()

    # Load MitoTip
    mitotip_map = {"Q1": "likely pathogenic", "Q2": "possibly pathogenic", "Q3": "possibly benign", "Q4": "likely benign"}
    mitotip_df = pd.read_csv(args.mitotipcache, sep="\t", dtype=str)
    mitotip_df["prediction"] = mitotip_df["Quartile"].map(mitotip_map)
    mitotip_db = pd.Series(
        mitotip_df["prediction"].values,
        index=list(zip(mitotip_df["rCRS"], mitotip_df["Position"], mitotip_df["Alt"]))
    ).to_dict()

    # Load HelixMTdb
    helix_df = pd.read_csv(args.helixcache, sep="\t", dtype=str)
    helix_df = helix_df[helix_df["alleles"].str.count(",") == 1].copy()
    helix_df["pos"] = helix_df["locus"].str.split("chrM:").str[1]
    a = helix_df["alleles"].str.split('"').str
    helix_df["ref"], helix_df["alt"] = a[1], a[3]
    for col in ["AF_hom", "max_ARF"]:
        helix_df[col] = pd.to_numeric(helix_df[col], errors="coerce").fillna(0.0)
    helix_df["max_het"] = helix_df.apply(lambda r: 1.0 if r["AF_hom"] > 0 else r["max_ARF"], axis=1)
    helix_db = pd.Series(
        list(zip(helix_df["max_het"], helix_df["AF_hom"], helix_df.get("AF_het", "0.0"))),
        index=pd.MultiIndex.from_frame(helix_df[["ref", "pos", "alt"]])
    ).to_dict()

    # -------------------- Step 3: Load and Parse VEP VCF ------------------- #
    with open(args.vep_vcf, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    vcf_df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str).rename(columns={"#CHROM": "CHROM"})
    
    sample_cols = list(vcf_df.columns[9:])
    if args.pipeline_mode == "disease" and proband_set:
        sample_cols = [c for c in sample_cols if any(c.startswith(p) for p in proband_set)]

    def extract_csq(info_str):
        if pd.isna(info_str): return ""
        for segment in str(info_str).split(";"):
            if segment.startswith("CSQ="): return segment[4:]
        return ""

    vcf_df["CSQ_STRING"] = vcf_df["INFO"].apply(extract_csq)
    vep_cols = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
        "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
        "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "HGVS_OFFSET",
    ]
    info_split = vcf_df["CSQ_STRING"].str.split("|", expand=True)
    for i, col in enumerate(vep_cols):
        vcf_df[col] = info_split.iloc[:, i] if i < info_split.shape[1] else pd.NA

    stacks = []
    for sc in sample_cols:
        sub = vcf_df.copy()
        sub["SampleID"] = sc
        sub["SAMPLE_DATA"] = sub[sc]
        stacks.append(sub)
    long_df = pd.concat(stacks, ignore_index=True)

    # -------------------- Step 4: Extract Genotype Fields ------------------ #
    fmt_headers = long_df["FORMAT"].iloc[0].split(":")
    sample_values = long_df["SAMPLE_DATA"].str.split(":", expand=True)

    def pick(tag: str):
        try:
            idx = fmt_headers.index(tag)
            return sample_values.iloc[:, idx]
        except (ValueError, IndexError):
            return pd.Series([pd.NA] * len(long_df))

    long_df["GT"] = pick("GT")
    long_df["DP"] = pd.to_numeric(pick("DP"), errors="coerce").fillna(0)
    long_df["AD"] = pick("AD")
    
    hl_raw = pick("HL")
    long_df["AF"] = pd.to_numeric(pd.Series(hl_raw), errors="coerce").fillna(0.0).values

    # -------------------- Step 5: Database Annotation --------------------- #
    var_keys = list(zip(long_df["REF"], long_df["POS"], long_df["ALT"]))
    
    # gnomAD: FIXED numeric conversion
    gvals = [gnomad_db.get(k, ("0", "0", "0")) for k in var_keys]
    long_df["gnomad_max_hl"] = [x[0] for x in gvals]
    long_df["gnomad_af_hom"] = pd.Series(pd.to_numeric([x[1] for x in gvals], errors="coerce")).fillna(0).values
    long_df["gnomad_af_het"] = pd.Series(pd.to_numeric([x[2] for x in gvals], errors="coerce")).fillna(0).values
    
    # Helix: FIXED numeric conversion
    hvals = [helix_db.get(k, (0.0, 0.0, 0.0)) for k in var_keys]
    long_df["helix_max_hl"] = [x[0] for x in hvals]
    long_df["helix_af_hom"] = pd.Series(pd.to_numeric([x[1] for x in hvals], errors="coerce")).fillna(0).values
    long_df["helix_af_het"] = pd.Series(pd.to_numeric([x[2] for x in hvals], errors="coerce")).fillna(0).values
    
    # ClinVar & MitoMap
    long_df["clinvar_interp"] = [clinvar_db.get(k, "") for k in var_keys]
    long_df["mitomap_gbcnt"]  = [mitomap_poly_db.get(k, ("0",))[0] for k in var_keys]
    
    # MitoMap AF
    long_df["mitomap_af"] = pd.Series(pd.to_numeric(long_df["mitomap_gbcnt"], errors="coerce")).fillna(0).values / 61134.0
    
    md = [mitomap_disease_db.get(k, ("", "", "", "")) for k in var_keys]
    long_df["mitomap_status"], long_df["mitomap_disease"] = [d[0] for d in md], [d[3] for d in md]
    long_df["mitomap_plasmy"] = [f"{d[1]}/{d[2]}" if d[1] or d[2] else "" for d in md]
    
    long_df["apogee_class"]   = [str(apogee_db.get(k, "")).strip("()").replace("'", "").replace(", ", "/") for k in var_keys]
    long_df["mitotip_class"]  = [mitotip_db.get(k, "") for k in var_keys]
    long_df["hmtvar_class"]   = [get_hmtvar(k, hmtvar_db) for k in var_keys]

    # NEW: nAPOGEE, tAPOGEE, MLC_score, and in_regional_constraint
    long_df["nAPOGEE"] = [napogee_db.get(k, ("",))[0] for k in var_keys]
    long_df["tAPOGEE"] = [tapogee_db.get(k, ("",))[0] for k in var_keys]

    def get_mlc_score(row):
        # 如果是 SNV，使用 REF+POS+ALT 匹配
        if row["VARIANT_CLASS"] == "SNV":
            key = (row["REF"], row["POS"], row["ALT"])
            return mlc_snv_db.get(key, ("",))[0]
        else:
            # 如果是 Indel，仅根据 POS 匹配
            return mlc_indel_db.get(row["POS"], "")

    long_df["MLC_score"] = long_df.apply(get_mlc_score, axis=1)

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

    long_df["in_regional_constraint"] = long_df.apply(get_rc_status, axis=1)

    # Stats
    long_df["variant_key"] = long_df["POS"].astype(str) + ":" + long_df["REF"] + ":" + long_df["ALT"]
    ac_counts = long_df[long_df["AF"] > 0]["variant_key"].value_counts().to_dict()
    long_df["in_cohort_AC"] = long_df["variant_key"].map(ac_counts).fillna(0).astype(int)
    long_df["Freq"] = (long_df["in_cohort_AC"] / len(sample_cols) * 100) if sample_cols else 0.0

    long_df = long_df[long_df["AF"] > 0].copy()

    # -------------------- Step 6: Export Results -------------------------- #
    base_cols = [
        "SampleID", "variant_key", "CHROM", "POS", "REF", "ALT", "FILTER", "GT", "AD", "AF", "DP", 
        "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "HGVSp", "Codons", "VARIANT_CLASS", 
        "gnomad_max_hl", "gnomad_af_hom", "gnomad_af_het", "apogee_class", "mitotip_class", 
        "hmtvar_class", "nAPOGEE", "tAPOGEE", "MLC_score", "in_regional_constraint",
        "helix_max_hl", "helix_af_hom", "helix_af_het", "mitomap_gbcnt", 
        "mitomap_af", "mitomap_status", "mitomap_plasmy", "mitomap_disease", "clinvar_interp", 
        "in_cohort_AC", "Freq"
    ]
    for c in base_cols:
        if c not in long_df.columns: long_df[c] = pd.NA
    
    pre_out = long_df[base_cols].rename(columns={"AF": "Heteroplasmy"}).copy()
    fname = "Proband_variant_list" if args.pipeline_mode == "disease" else "variant_list"
    pre_out.to_csv(os.path.join(outdir, f"{fname}_prefiltering.txt"), sep="\t", index=False, na_rep="")

    for col in ["gnomad_af_hom", "helix_af_hom", "mitomap_af", "Heteroplasmy", "Freq"]:
        pre_out[col] = pd.to_numeric(pre_out[col], errors="coerce")

    filt = pre_out[
        (pre_out["gnomad_af_hom"] < 0.01) &
        (pre_out["helix_af_hom"] < 0.01) &
        (pre_out["mitomap_af"] < 0.01) &
        (pre_out["Consequence"] != "synonymous_variant") &
        (pre_out["Freq"] < 100) &
        (pre_out["Heteroplasmy"] > 0.01) & 
        (~pre_out["clinvar_interp"].isin(["Benign", "Likely benign"]))
    ].copy()

    if not filt.empty:
        filt["__has_disease__"] = filt["mitomap_disease"].apply(lambda x: 0 if pd.isna(x) or x == "" else 1)
        filt = filt.sort_values(by=["__has_disease__", "gnomad_af_hom", "helix_af_hom", "mitomap_af"], ascending=[False, True, True, True]).drop(columns="__has_disease__")

    filt.to_csv(os.path.join(outdir, f"{fname}.txt"), sep="\t", index=False, na_rep="")
    print(f"[+] Final tables written to {outdir}. Total records found: {len(long_df)}")

if __name__ == "__main__":
    main()