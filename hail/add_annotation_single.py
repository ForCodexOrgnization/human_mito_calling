#!/usr/bin/env python3

from logging import config
import os
import sys
import json
import csv
import glob
import subprocess
import pandas as pd

# ==============================================================================
# Helpers
# ==============================================================================

def run_command(command, cwd=None, is_shell=False):

    cmd_str = command if is_shell else " ".join(command)
    print(f"[*] Executing (in {cwd or '.'}): {cmd_str}")
    process = subprocess.Popen(
        cmd_str if is_shell else command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=cwd,
        shell=is_shell,
        executable="/bin/bash" if is_shell else None,
    )
    for line in process.stdout:
        print(line, end="")
    rc = process.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd_str)
    print(f"[+] Command successful: {cmd_str}")


def read_files_to_dataframe(pattern, extractor, col_names):

    frames = []
    for filepath in glob.glob(pattern):
        with open(filepath, "r") as f:
            df = pd.read_csv(f, sep="\t", dtype=str)
        frames.append(extractor(df))
    if not frames:
        return pd.DataFrame(columns=col_names)
    combined = pd.concat(frames, ignore_index=True)
    combined.columns = col_names
    return combined


def _run_vep_single(file_path, config):
    """
    Run VEP once for a single VCF file using config-specified cache/output dirs.
    """
    try:
        out_dir = config["vep_vcf_dir"]
        cache_dir = config["vep_cache_dir"]
        os.makedirs(out_dir, exist_ok=True)

        base = os.path.basename(file_path).replace(".final.split.vcf", "")
        out_vcf = os.path.join(out_dir, f"{base}_vep.vcf")

        print(f"[*] Starting VEP for: {base}")
        vep_args = [
            "vep", "--cache",
            "-i", file_path,
            "-o", out_vcf,
            "--dir_cache", cache_dir,
            "--no_stats",
            "--distance", "0",
            "--biotype",
            "--symbol",
            "--hgvs",
            "--variant_class",
            "--force_overwrite",
            "--vcf",
        ]
        run_command(vep_args, is_shell=False)
        print(f"[+] Finished VEP for: {base}")
    except Exception as e:
        raise RuntimeError(f"VEP failed for {file_path}: {e}")

# ==============================================================================
# Main single-sample steps
# ==============================================================================

def step1_run_vep_single(config):
    """
    Pick the input VCF and run VEP.

    Priority:
      - config['single_vcf_path'] if present and exists
      - otherwise, the first *.final.split.vcf under config['raw_vcf_dir']
    """
    raw_dir = config["raw_vcf_dir"]
    os.makedirs(config["vep_vcf_dir"], exist_ok=True)

    vcf_path = config.get("single_vcf_path")
    if vcf_path and os.path.isfile(vcf_path):
        target = vcf_path
    else:
        matches = sorted(glob.glob(os.path.join(raw_dir, "*.final.split.vcf")))
        if not matches:
            print(f"[!] No VCF found in {raw_dir}. Skipping VEP.")
            return
        target = matches[0]

    print(f"[*] Using VCF: {target}")
    _run_vep_single(target, config)


def step2_extract_metadata(config):
    """
    Extract haplogroup/contamination and write standardized metadata TSVs.
    """
    print("\n--- Step 2: Extract haplogroup and contamination ---")
    source_dir = config["haplocheck_dir"]
    out_dir = config["metadata_dir"]
    os.makedirs(out_dir, exist_ok=True)

    file_pattern = os.path.join(source_dir, "*.txt")

    # Haplogroup
    haplogroup_df = read_files_to_dataframe(
        file_pattern,
        lambda df: df[["SampleID", "HgMajor"]],
        ["Sample_ID", "Full_Haplogroup"],
    )
    hap_path = os.path.join(out_dir, "haplogroup_full.txt")
    haplogroup_df.to_csv(hap_path, sep="\t", index=False)
    print(f"[+] Haplogroup file saved to: {hap_path}")

    # Contamination
    contam_df = read_files_to_dataframe(
        file_pattern,
        lambda df: df[["SampleID", "Contamination"]],
        ["Sample_ID", "Contamination_Status"],
    )
    contam_path = os.path.join(out_dir, "contamination.txt")
    contam_df.to_csv(contam_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    print(f"[+] Contamination file saved to: {contam_path}")


def step3_run_variant_processor(config):
    """
    Invoke the single-sample variant post-processing script.
    """
    print("\n--- Step 3: Run single-sample variant processing ---")
    cmd = [
        sys.executable, "./hail/process_variants_single.py",
        "--vep-vcf-dir",        config["vep_vcf_dir"],
        "--final-output-dir",   config["final_output_dir"],
        "--fullhaplogroups",    config["fullhaplogroups"],
        "--contamination",      config["contamination"],
        "--gnomadcache",        config["gnomadcache"],
        "--clinvarcache",       config["clinvarcache"],
        "--mitomap-polycache",  config["mitomap_polycache"],
        "--mitomap-diseasecache", config["mitomap_diseasecache"],
        "--helixcache",         config["helixcache"],
        "--haplogroup-varcache", config["haplogroup_varcache"],
        "--mitimpactcache",     config["mitimpactcache"],
        "--mitotipcache",       config["mitotipcache"],
        "--hmtvarcache",        config["hmtvarcache"],
        # ---- NEW annotation caches ----
        "--napogeecache",       config["napogeecache"],
        "--tapogeecache",       config["tapogeecache"],
        "--mlc_snv_cache",      config["mlc_snv_cache"],
        "--mlc_indel_cache",    config["mlc_indel_cache"],
        "--regional_constraint_cache", config["regional_constraint_cache"],
    ]

    # ---- OPTIONAL disease meta ----
    disease_meta = config.get("disease_meta_file", None)
    if disease_meta is not None:
        disease_meta = str(disease_meta).strip()
        if disease_meta and disease_meta.lower() not in ("none", "null"):
            if os.path.isfile(disease_meta):
                cmd += ["--disease-meta-file", disease_meta]
            else:
                print(f"[WARN] disease_meta_file not found, skip: {disease_meta}")

    run_command(cmd, is_shell=False)


# ==============================================================================
# Entrypoint
# ==============================================================================

def main():
    with open("config.json", "r") as f:
        config = json.load(f)

    # Ensure output structure exists
    os.makedirs(config["vep_vcf_dir"], exist_ok=True)
    os.makedirs(config["metadata_dir"], exist_ok=True)
    os.makedirs(config["final_output_dir"], exist_ok=True)

    try:
        step1_run_vep_single(config)
        step2_extract_metadata(config)
        step3_run_variant_processor(config)
        print("\n[SUCCESS] Single-sample pipeline completed.")
    except Exception as e:
        print(f"\n[ERROR] Pipeline terminated: {e}")
        import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()