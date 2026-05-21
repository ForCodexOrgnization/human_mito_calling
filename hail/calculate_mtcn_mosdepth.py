#!/usr/bin/env python3
"""
calculate_mtcn_v3.py
====================
Calculate Mitochondrial DNA Copy Number (mtCN) using mosdepth.
Formula: mtCN_final = 2 × mean(mtDNA_coverage) / median(nuclear_coverage)
Note: mosdepth path is passed as a command-line argument.
"""

import os
import argparse
import subprocess
import pandas as pd
import numpy as np

# -------------------------------------------------------------------------
# Step 1 & 2: Generate and parse nuclear coverage using mosdepth
# -------------------------------------------------------------------------
def get_nuclear_coverage_mosdepth(cram, ref_fasta, output_dir, sample_id, mosdepth_path):
    """
    Run mosdepth to obtain whole-genome coverage statistics.
    Returns: (nuc_median, nuc_mean)
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, sample_id)
    
    # Construct mosdepth command
    # -n: Do not output per-base files (extremely fast)
    # -Q 20: Filter reads with Mapping Quality < 20
    # --fast-mode: Enable further acceleration
    cmd = [
        mosdepth_path, "-n", 
        "-Q", "20", 
        "--fast-mode",
        "-f", ref_fasta,
        prefix, cram
    ]

    print(f"[CMD] {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERR] mosdepth execution failed: {e}")
        raise RuntimeError(f"mosdepth failed for {sample_id}")

    # Parse summary file: <prefix>.mosdepth.summary.txt
    summary_file = f"{prefix}.mosdepth.summary.txt"
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")

    # Parse whole-genome mean (from the 'total' row)
    summary_df = pd.read_csv(summary_file, sep="\t")
    nuc_mean = float(summary_df.loc[summary_df['chrom'] == 'total', 'mean'].values[0])

    # Parse median from <prefix>.mosdepth.global.dist.txt (at cumulative fraction 0.50)
    dist_file = f"{prefix}.mosdepth.global.dist.txt"
    nuc_median = 0.0
    with open(dist_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == "total":
                cum_frac = float(parts[2])
                if cum_frac >= 0.50:
                    nuc_median = float(parts[1])
                    break
    
    print(f"[OK] Nuclear (mosdepth): mean={nuc_mean:.2f}, median={nuc_median:.2f}")
    return nuc_median, nuc_mean


# -------------------------------------------------------------------------
# Step 3: Calculate mitochondrial coverage from per_base_coverage.tsv
# -------------------------------------------------------------------------
def calculate_mt_coverage(mt_coverage_file):
    """Calculate mean and median mtDNA coverage from per-base TSV"""
    if not os.path.exists(mt_coverage_file):
        raise FileNotFoundError(f"Mitochondrial coverage file not found: {mt_coverage_file}")

    df = pd.read_csv(mt_coverage_file, sep="\t")
    # Automatically identify column name (usually 'coverage' or 'depth')
    cov_col = "coverage" if "coverage" in df.columns else "depth"
    
    cov_mean = float(pd.to_numeric(df[cov_col], errors="coerce").mean())
    cov_median = float(pd.to_numeric(df[cov_col], errors="coerce").median())
    
    print(f"[OK] mtDNA: mean={cov_mean:.2f}, median={cov_median:.2f}")
    return cov_mean, cov_median


# -------------------------------------------------------------------------
# Step 4: Main Control Logic
# -------------------------------------------------------------------------
def main(args):
    sample_id = os.path.basename(args.cram).replace(".cram", "")
    os.makedirs(args.output, exist_ok=True)

    # 1. Get nuclear genome depth using the provided mosdepth path
    nuc_median, nuc_mean = get_nuclear_coverage_mosdepth(
        args.cram, args.ref_fasta, args.output, sample_id, args.mosdepth
    )

    # 2. Get mitochondrial genome depth
    mt_mean, mt_median = calculate_mt_coverage(args.mt_coverage)

    # 3. Compute mtCN metrics
    # Formula: 2 * Mean Mitochondrial Coverage / Median Nuclear Coverage
    mtcn_final = 2 * mt_mean / nuc_median if (nuc_median > 0) else np.nan

    contam = 0.0  # Placeholder for contamination estimation (if implemented later)
    is_contaminated_str = "NO"

    # Apply QC filtering logic

    if hasattr(args, 'haplocheck') and args.haplocheck and os.path.exists(args.haplocheck):
        try:
            h_df = pd.read_csv(args.haplocheck, sep="\t")
            if 'Contamination' in h_df.columns:
                is_contaminated_str = str(h_df['Contamination'].values[0]).upper()
                
                if is_contaminated_str == "YES":
                    contam = 1.0  
                else:
                    contam = 0.0
        except Exception as e:
            print(f"[WARN] Failed to parse haplocheck file: {e}")

    pass_filter = (
        (mt_mean >= args.min_mt_cov)
        and (np.isfinite(mtcn_final))
        and (args.min_mtcn <= mtcn_final <= args.max_mtcn)
        and (contam <= args.max_contam)
        and (is_contaminated_str != "YES")
    )

    # 4. Summarize records
    record = {
        "Sample_ID": sample_id,
        "Mean_mtDNA_Coverage": round(mt_mean, 2),
        "Median_mtDNA_Coverage": round(mt_median, 2),
        "Mean_Nuclear_Coverage": round(nuc_mean, 2),
        "Median_Nuclear_Coverage": round(nuc_median, 2),
        "mtCN_final": round(mtcn_final, 2),
        "Contamination": is_contaminated_str,
        "Pass_Filter": pass_filter
    }

    out_tsv = os.path.join(args.output, "mtCN_summary.txt")
    pd.DataFrame([record]).to_csv(out_tsv, sep="\t", index=False)
    
    print(f"[DONE] {sample_id} processed successfully. mtCN: {mtcn_final:.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mtCN Calculation (mosdepth version)")
    # Required parameters
    parser.add_argument("--cram", required=True, help="Input CRAM file path")
    parser.add_argument("--crai", required=True, help="Index CRAI file (not used here; kept for interface parity).")
    parser.add_argument("--ref_fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--mt_coverage", required=True, help="Mitochondrial per-base coverage TSV")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--mosdepth", required=True, help="Path to mosdepth executable")
    parser.add_argument("--haplocheck", help="Haplocheck contamination file (optional)")
    
    # Optional QC parameters
    parser.add_argument("--min_mt_cov", type=int, default=100, help="Minimum mtDNA coverage threshold")
    parser.add_argument("--min_mtcn", type=float, default=50, help="Minimum mtCN threshold")
    parser.add_argument("--max_mtcn", type=float, default=500, help="Maximum mtCN threshold")
    parser.add_argument("--max_contam", type=float, default=0.02, help="Maximum contamination allowed.")


    args = parser.parse_args()
    main(args)