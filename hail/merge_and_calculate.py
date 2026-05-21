#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Merge per-sample results and calculate cohort AC/AF")
    parser.add_argument("--input-dir", required=True, help="Directory with sample.pre.txt files")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--pipeline-mode", default="population", help="Mode: population or disease")
    parser.add_argument("--do_post_filter", type=str, default="True", help="Whether to apply strict filters (True/False)")
    args = parser.parse_args()

    do_post_filter = args.do_post_filter.lower() == "true"

    all_dfs = []
    files = [f for f in os.listdir(args.input_dir) if f.endswith(".pre.txt")]
    
    if not files:
        print(f"[!] No .pre.txt files found in {args.input_dir}. Check your input directory.")
        return

    print(f"[*] Merging {len(files)} sample files...")
    for f in files:
        df = pd.read_csv(os.path.join(args.input_dir, f), sep="\t", dtype=str)
        all_dfs.append(df)
    
    combined_df = pd.concat(all_dfs, ignore_index=True)

    num_samples = len(files)

    print("[*] Calculating cohort-level AC and AF...")
    combined_df["variant_key"] = (
        combined_df["POS"].astype(str) + ":" + 
        combined_df["REF"].astype(str) + ":" + 
        combined_df["ALT"].astype(str)
    )
    
    ac_map = combined_df["variant_key"].value_counts().to_dict()
    combined_df["in_cohort_AC"] = combined_df["variant_key"].map(ac_map)

    print("[*] Calculating common_low_heteroplasmy flag...")
    temp_hl = pd.to_numeric(combined_df["Heteroplasmy"], errors="coerce").fillna(0.0)
    
    mid_het_mask = (temp_hl > 0) & (temp_hl < 0.5) & (combined_df["FILTER"].isin(["PASS", ".", ""]))
    
    ac_mid_het_map = combined_df[mid_het_mask]["variant_key"].value_counts().to_dict()
    
    combined_df["common_low_heteroplasmy"] = combined_df["variant_key"].apply(
        lambda x: (ac_mid_het_map.get(x, 0) / num_samples) > 0.001
    )

    combined_df["in_cohort_AF"] = combined_df["in_cohort_AC"].astype(float) / num_samples

    prefix = "Proband_" if args.pipeline_mode == "disease" else ""
    os.makedirs(args.out_dir, exist_ok=True)
    pre_out = os.path.join(args.out_dir, f"{prefix}variant_list_prefiltering.txt")
    final_out = os.path.join(args.out_dir, f"{prefix}variant_list.txt")

    combined_df.drop(columns=["variant_key"], inplace=True)
    combined_df.to_csv(pre_out, sep="\t", index=False)
    print(f"[+] Merged prefiltering table saved to: {pre_out}")

    if do_post_filter:
        print("[*] do_post_filter is True: Applying strict cohort filters...")
        
        numeric_cols = ["gnomad_af_hom", "helix_af_hom", "mitomap_af", "Heteroplasmy", "in_cohort_AF"]
        for col in numeric_cols:
            if col in combined_df.columns:
                combined_df[col] = pd.to_numeric(combined_df[col], errors="coerce").fillna(0.0)
        
        mask = (
            (~combined_df["FILTER"].str.contains("indel_stack", na=False)) & 
            (combined_df["FILTER"].isin(["PASS", ".", ""])) &
            (combined_df["Haplogroup_Var_Status"] != "haplo_var_match")
        )
        
        if "Consequence" in combined_df.columns:
            mask &= (~combined_df["Consequence"].str.contains("synonymous", na=False))

        filtered_df = combined_df[mask].copy()
        filtered_df.to_csv(final_out, sep="\t", index=False)
        print(f"[+] Filtered variant list saved to: {final_out}")

    print(f"[SUCCESS] Cohort processing finished.")

if __name__ == "__main__":
    main()