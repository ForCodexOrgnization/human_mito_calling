#!/usr/bin/env python3

import os
import io
import json
import argparse
import pandas as pd
import numpy as np
import gc
import re

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def load_complex_db(path, key_cols, val_cols, delimiter="\t"):
    if not os.path.exists(path):
        print(f"WARNING: Database {path} not found.")
        return {}
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = df.columns.str.strip()
    idx = pd.MultiIndex.from_frame(df[key_cols])
    df["__values__"] = df[val_cols].apply(tuple, axis=1)
    return pd.Series(df["__values__"].values, index=idx).to_dict()

def load_pos_only_db(path, pos_col, val_col, delimiter="\t"):
    if not os.path.exists(path): return {}
    df = pd.read_csv(path, sep=delimiter, low_memory=False, dtype=str, comment="#", encoding='ISO-8859-1')
    df.columns = df.columns.str.strip()
    return pd.Series(df[val_col].values, index=df[pos_col]).to_dict()

def get_hmtvar(k, hmtvar_db):
    val = hmtvar_db.get(k)
    if val and val[0]:
        try: return json.loads(val[0]).get("pathogenicity", "")
        except: return ""
    return ""

def load_haplo_assoc(path):
    if not os.path.exists(path): 
        print(f"WARNING: HaploVar database {path} not found.")
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    return pd.Series(df.iloc[:, 1].str.lower().values, index=df.iloc[:, 0]).to_dict()

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep-vcf", required=True)
    ap.add_argument("--final-output-dir", required=True)
    ap.add_argument("--gnomadcache", required=True)
    ap.add_argument("--clinvarcache", required=True)
    ap.add_argument("--mitomap-polycache", required=True)
    ap.add_argument("--mitomap-diseasecache", required=True)
    ap.add_argument("--helixcache", required=True)
    ap.add_argument("--mitimpactcache", required=True)
    ap.add_argument("--mitotipcache", required=True)
    ap.add_argument("--hmtvarcache", required=True)
    ap.add_argument("--napogeecache", required=True)
    ap.add_argument("--tapogeecache", required=True)
    ap.add_argument("--mlc_snv_cache", required=True)
    ap.add_argument("--mlc_indel_cache", required=True)
    ap.add_argument("--regional_constraint_cache", required=True)
    ap.add_argument("--pipeline-mode", choices=["population", "disease"], required=True)
    ap.add_argument("--fullhaplogroups", required=True)
    ap.add_argument("--haplovarcache", required=True)
    ap.add_argument("--disease-meta-file")

    args = ap.parse_args()
    outdir = args.final_output_dir
    os.makedirs(outdir, exist_ok=True)

    print("[*] Loading Databases...")
    gnomad_db = load_complex_db(args.gnomadcache, ["ref", "position", "alt"], ["max_observed_heteroplasmy", "AF_hom", "AF_het"])
    mitomap_poly_db = load_complex_db(args.mitomap_polycache, ["ref", "pos", "alt"], ["gbcnt"])
    mitomap_disease_db = load_complex_db(args.mitomap_diseasecache, ["ref", "pos", "alt"], ["status", "homoplasmy", "heteroplasmy", "disease"])
    apogee_db = load_complex_db(args.mitimpactcache, ["Ref", "Start", "Alt"], ["APOGEE1", "APOGEE2"])
    hmtvar_db = load_complex_db(args.hmtvarcache, ["REF", "POS", "ALT"], ["HmtVar"])
    haplovar_db = load_haplo_assoc(args.haplovarcache)
    napogee_db = load_complex_db(args.napogeecache, ["ref", "start", "alt"], ["nAPOGEE_score"])
    tapogee_db = load_complex_db(args.tapogeecache, ["Ref", "Pos", "Alt"], ["t-APOGEE score"])
    mlc_snv_db = load_complex_db(args.mlc_snv_cache, ["Reference", "Position", "Alternate"], ["MLC_score"])
    mlc_indel_db = load_pos_only_db(args.mlc_indel_cache, "Position", "MLC_pos_score")

    haplo_map = {}
    if args.fullhaplogroups and os.path.isdir(args.fullhaplogroups):
        print(f"[*] Reading haplocheck files from directory: {args.fullhaplogroups}")
        for f_name in os.listdir(args.fullhaplogroups):
            if f_name.endswith(".haplocheck_contamination.txt"):
                sample_id_from_file = f_name.replace(".haplocheck_contamination.txt", "")
                
                file_path = os.path.join(args.fullhaplogroups, f_name)
                try:
                    with open(file_path, "r") as f:
                        valid_lines = [l.strip() for l in f if l.strip()]
                        if len(valid_lines) > 1:
                            header = valid_lines[0].split("\t")
                            data = valid_lines[1].split("\t")
        
                        if "HgMajor" in header:
                            hg_idx = header.index("HgMajor")
                            haplo_map[sample_id_from_file] = data[hg_idx]
                            print(f"DEBUG: Loaded {sample_id_from_file} -> {data[hg_idx]}")

                except Exception as e:
                    print(f"WARNING: Could not read {f_name}: {e}")
    
    elif args.fullhaplogroups and os.path.isfile(args.fullhaplogroups):
        with open(args.fullhaplogroups, "r") as f:
            for line in f:
                if line.startswith("Sample") or not line.strip(): continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    haplo_map[parts[0]] = parts[1]
    
    rc_df = pd.read_csv(args.regional_constraint_cache, sep="\t", dtype=str, encoding='ISO-8859-1')
    rc_db = pd.Series(list(zip(rc_df["in_rc"], rc_df["min_distance_to_rc"])), index=pd.MultiIndex.from_frame(rc_df[["REF", "POS", "ALT"]])).to_dict()
    del rc_df

    cv_df = pd.read_csv(args.clinvarcache, sep="\t", dtype=str, encoding='ISO-8859-1')
    cv_df.columns = cv_df.columns.str.strip()
    spdi = cv_df["Canonical SPDI"].str.split(":", expand=True)
    clinvar_db = pd.Series(cv_df["Germline classification"].values, index=pd.MultiIndex.from_frame(pd.concat([spdi[2], cv_df["GRCh38Location"], spdi[3]], axis=1))).to_dict()
    del cv_df

    mtip_df = pd.read_csv(args.mitotipcache, sep="\t", dtype=str)
    mtip_map = {"Q1": "pathogenic", "Q2": "possibly pathogenic", "Q3": "possibly benign", "Q4": "benign"}
    mitotip_db = pd.Series(mtip_df["Quartile"].map(mtip_map).values, index=list(zip(mtip_df["rCRS"], mtip_df["Position"], mtip_df["Alt"]))).to_dict()
    del mtip_df

    hx_df = pd.read_csv(args.helixcache, sep="\t", dtype=str)
    hx_df.columns = hx_df.columns.str.strip()
    hx_df["pos"] = hx_df["locus"].str.extract(r'chrM:(\d+)')
    all_bases = hx_df["alleles"].str.findall(r'([A-Za-z]+)')
    hx_df = hx_df[all_bases.str.len() == 2].copy()
    
    hx_df["ref"] = all_bases.str[0]
    hx_df["alt"] = all_bases.str[1]

    for c in ["AF_hom", "AF_het", "max_ARF"]:
        hx_df[c] = pd.to_numeric(hx_df.get(c, "0.0"), errors="coerce").fillna(0.0)

    hx_df["max_het"] = hx_df.apply(
        lambda r: 1.0 if r["AF_hom"] > 0 else r["max_ARF"], axis=1
    )

    helix_db = pd.Series(
        list(zip(hx_df["max_het"], hx_df["AF_hom"], hx_df["AF_het"])), 
        index=pd.MultiIndex.from_frame(hx_df[["ref", "pos", "alt"]])
    ).to_dict()
    del hx_df; gc.collect()

    with open(args.vep_vcf, "r") as f:

        lines = [l for l in f if not l.startswith("##")]

    vcf_df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    vcf_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
    
    vcf_df["CSQ_STR"] = vcf_df["INFO"].str.extract(r'CSQ=([^;]+)')
    vep_cols = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_pos", "CDS_pos", "Prot_pos", "Amino_acids", "Codons", "Existing_var", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS"]
    csq_split = vcf_df["CSQ_STR"].str.split("|", expand=True)
    for i, col in enumerate(vep_cols):
        vcf_df[col] = csq_split.iloc[:, i] if i < csq_split.shape[1] else pd.NA
    vcf_df.drop(columns=["INFO", "CSQ_STR"], inplace=True)

    # 6. 流式处理各样本
    sample_cols = list(vcf_df.columns[9:])
    fmt_headers = vcf_df["FORMAT"].iloc[0].split(":")
    
    af_tag = next((t for t in ["HL", "AF"] if t in fmt_headers), None)
    af_idx = fmt_headers.index(af_tag) if af_tag else -1
    gt_idx = fmt_headers.index("GT") if "GT" in fmt_headers else -1
    ad_idx = fmt_headers.index("AD") if "AD" in fmt_headers else -1
    dp_idx = fmt_headers.index("DP") if "DP" in fmt_headers else -1

    # --- 新增/优化: Indel Stack 统计逻辑 ---
    print("[*] Calculating indel_stack metrics (optimized for split VCF)...")
    pos_sample_indels = {}
    
    def safe_get_vector(idx, df_len, sample_vals_df, default_val):
        if idx != -1 and idx < sample_vals_df.shape[1]:
            return sample_vals_df.iloc[:, idx]
        return pd.Series([default_val] * df_len)

    for sc in sample_cols:
        s_vals = vcf_df[sc].str.split(":", expand=True)
        s_af = pd.to_numeric(safe_get_vector(af_idx, len(vcf_df), s_vals, "0.0"), errors="coerce").fillna(0.0)
        
        # 判定 Indel: 只要 REF 长度 != ALT 长度，就是 Indel (比 VARIANT_CLASS 更准)
        is_indel = vcf_df["REF"].str.strip().str.len() != vcf_df["ALT"].str.strip().str.len()
        mask = (s_af > 0.01) & is_indel
        if not mask.any(): continue
        
        # 记录样本在该 POS 下拥有的所有 ALT
        for _, row in vcf_df[mask].iterrows():
            p = str(row["POS"]).strip()
            ref = str(row["REF"]).strip()
            # 即使是拆分后的 VCF，防御性处理逗号
            alts = str(row["ALT"]).strip().split(",")
            for a in alts:
                a = a.strip()
                if len(ref) != len(a):
                    if p not in pos_sample_indels: pos_sample_indels[p] = {}
                    if sc not in pos_sample_indels[p]: pos_sample_indels[p][sc] = set()
                    pos_sample_indels[p][sc].add(a)

    # 判定变异级别的 stack 状态
    variant_stack_status = {}
    for pos, s_map in pos_sample_indels.items():
        all_alts_at_pos = set().union(*s_map.values())
        for alt in all_alts_at_pos:
            # 找到所有携带该 (pos, alt) 的样本
            carrying_samples = [s for s, sample_alts in s_map.items() if alt in sample_alts]
            # 标记规则：如果在所有携带该变异的样本中，该位置都存在 >= 2 个不同的 Indel
            if carrying_samples and all(len(s_map[s]) >= 2 for s in carrying_samples):
                variant_stack_status[(pos, alt)] = True
    # ---------------------------------

    fname = "Proband_variant_list" if args.pipeline_mode == "disease" else "variant_list"
    pre_out_path = os.path.join(outdir, f"{fname}_prefiltering.txt")
    base_cols = ["SampleID", "variant_key", "CHROM", "POS", "REF", "ALT", "FILTER", "GT", "AD", "Heteroplasmy", "DP", "Consequence", "SYMBOL", "BIOTYPE", "HGVSc", "HGVSp", "VARIANT_CLASS", "indel_stack_filter","Haplogroup", "Haplogroup_Var_Status","gnomad_max_hl", "gnomad_af_hom", "gnomad_af_het", "apogee_class", "mitotip_class", "hmtvar_class", "nAPOGEE", "tAPOGEE", "MLC_score", "in_regional_constraint", "helix_max_hl", "helix_af_hom", "helix_af_het", "mitomap_gbcnt", "mitomap_af", "mitomap_status", "mitomap_plasmy", "mitomap_disease", "clinvar_interp"]

    first_write = True
    print(f"[*] Processing {len(sample_cols)} samples with robustness checks...")

    for sc in sample_cols:
        sample_vals = vcf_df[sc].str.split(":", expand=True)
        num_fields = sample_vals.shape[1]
        
        def safe_get(idx, default_val):
            if idx != -1 and idx < num_fields:
                return sample_vals.iloc[:, idx]
            return pd.Series([default_val] * len(vcf_df))

        cur_af_raw = safe_get(af_idx, "0.0")
        cur_af = pd.to_numeric(cur_af_raw, errors="coerce").fillna(0.0)
        
        mask = cur_af > 0.01
        if not mask.any(): continue
        
        sub = vcf_df[mask].copy().reset_index(drop=True)

        sub["SampleID"] = sc
        
        sub["FILTER"] = vcf_df.loc[mask, "FILTER"].values

        sub["Haplogroup"] = sub["SampleID"].map(haplo_map).fillna("")

        def get_haplo_status(row):
            key = f"{row['REF']}{row['POS']}{row['ALT']}"
            assoc = haplovar_db.get(key)
            if assoc:
                samp_h = str(row.get("Haplogroup", "")).lower()
                return "haplo_var_match" if samp_h in assoc else "haplo_var_diff_haplo"
            return "not_haplo_var"

        sub["Haplogroup_Var_Status"] = sub.apply(get_haplo_status, axis=1)

        sub["Heteroplasmy"] = cur_af[mask].values
        sub["GT"] = safe_get(gt_idx, "0/1")[mask].values
        sub["AD"] = safe_get(ad_idx, ".")[mask].values
        sub["DP"] = safe_get(dp_idx, ".")[mask].values
        sub["variant_key"] = sub["POS"] + ":" + sub["REF"] + ":" + sub["ALT"]
        
        sub["indel_stack_filter"] = sub.apply(lambda r: "indel_stack" if variant_stack_status.get((str(r["POS"]).strip(), str(r["ALT"]).strip()), False) else "pass", axis=1)

        keys = list(zip(sub["REF"], sub["POS"], sub["ALT"]))
        
        g = [gnomad_db.get(k, ("0","0","0")) for k in keys]
        sub["gnomad_max_hl"], sub["gnomad_af_hom"], sub["gnomad_af_het"] = [x[0] for x in g], [x[1] for x in g], [x[2] for x in g]
        
        h = [helix_db.get(k, ("0.0","0.0","0.0")) for k in keys]
        sub["helix_max_hl"], sub["helix_af_hom"], sub["helix_af_het"] = [x[0] for x in h], [x[1] for x in h], [x[2] for x in h]
        
        sub["clinvar_interp"] = [clinvar_db.get(k, "") for k in keys]
        sub["mitomap_gbcnt"] = [mitomap_poly_db.get(k, ("0",))[0] for k in keys]
        sub["mitomap_af"] = pd.to_numeric(sub["mitomap_gbcnt"], errors="coerce").fillna(0) / 61134.0
        
        md = [mitomap_disease_db.get(k, ("", "", "", "")) for k in keys]
        sub["mitomap_status"], sub["mitomap_disease"] = [d[0] for d in md], [d[3] for d in md]
        sub["mitomap_plasmy"] = [f"{d[1]}/{d[2]}" if d[1] or d[2] else "" for d in md]
        
        sub["apogee_class"] = [str(apogee_db.get(k, "")).strip("()").replace("'", "").replace(", ", "/") for k in keys]
        sub["mitotip_class"] = [mitotip_db.get(k, "") for k in keys]
        sub["hmtvar_class"] = [get_hmtvar(k, hmtvar_db) for k in keys]
        sub["nAPOGEE"] = [napogee_db.get(k, ("",))[0] for k in keys]
        sub["tAPOGEE"] = [tapogee_db.get(k, ("",))[0] for k in keys]
        sub["MLC_score"] = sub.apply(lambda r: mlc_snv_db.get((r["REF"], r["POS"], r["ALT"]), ("",))[0] if r["VARIANT_CLASS"] == "SNV" else mlc_indel_db.get(r["POS"], ""), axis=1)

        def get_rc(row):
            if "missense" not in str(row["Consequence"]) and "rrna" not in str(row["Consequence"]): return "NA"
            info = rc_db.get((row["REF"], row["POS"], row["ALT"]))
            if not info: return "no"
            return "yes" if str(info[0]).lower() == "yes" else ("proximal" if float(info[1]) <= 6 else "no")
        sub["in_regional_constraint"] = sub.apply(get_rc, axis=1)

        sub[base_cols].to_csv(pre_out_path, sep="\t", index=False, mode=('w' if first_write else 'a'), header=first_write)
        first_write = False
        del sub, sample_vals; gc.collect()

    if not os.path.exists(pre_out_path):
        print("[!] No records found."); return

    final_df = pd.read_csv(pre_out_path, sep="\t", low_memory=False)
    counts = final_df["variant_key"].value_counts().to_dict()
    final_df["in_cohort_AC"] = final_df["variant_key"].map(counts)
    final_df["in_cohort_AF"] = final_df["in_cohort_AC"] / len(sample_cols)

    num_cols = ["gnomad_af_hom", "helix_af_hom", "mitomap_af", "Heteroplasmy", "in_cohort_AF"]
    for col in num_cols:
        final_df[col] = pd.to_numeric(final_df[col], errors="coerce").fillna(0.0)

    clinvar_junk_pos = ["73", "263", "16159", "16182", "16183", "16223"]

    filt = final_df[
        (final_df["gnomad_af_hom"] < 0.01) &
        (final_df["helix_af_hom"] < 0.01) &
        (final_df["mitomap_af"] < 0.01) &
        (final_df["Consequence"].str.contains("synonymous") == False) &
        (final_df["in_cohort_AF"] <= 1.0) &
        (final_df["Heteroplasmy"] > 0.01) & 
        (final_df["indel_stack_filter"] != "indel_stack") & 
        (final_df["FILTER"].isin(["PASS", "."])) &
        (~final_df["clinvar_interp"].fillna("").isin(["Benign", "Likely benign", "Affects"])) &
        (~final_df["POS"].astype(str).isin(clinvar_junk_pos))
    ].copy()

    filt["__patho__"] = filt["mitomap_disease"].apply(lambda x: 0 if pd.isna(x) or x == "" else 1)
    filt = filt.sort_values(by=["__patho__", "gnomad_af_hom"], ascending=[False, True]).drop(columns="__patho__")

    filt.to_csv(os.path.join(outdir, f"{fname}.txt"), sep="\t", index=False)
    print(f"[SUCCESS] Processed {len(sample_cols)} samples. Output in {outdir}")

if __name__ == "__main__":
    main()