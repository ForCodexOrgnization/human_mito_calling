#!/usr/bin/env python3
import sys, json, argparse, subprocess, csv
from pathlib import Path

# -------------------- Utilities --------------------
def run(cmd, cwd=None):
    """Run a subprocess, stream stdout, raise on non-zero exit."""
    print(f"[*] Executing: {' '.join(map(str, cmd))}")
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True)
    for line in p.stdout:
        print(line, end="")
    rc = p.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)

# Required configuration keys (Updated with 5 new keys)
REQ_KEYS = [
    "vep_cache_dir", "pipeline_mode",
    "gnomadcache", "clinvarcache",
    "mitomap_polycache", "mitomap_diseasecache",
    "helixcache",
    "mitimpactcache", "mitotipcache", "hmtvarcache",
    "napogeecache", "tapogeecache", "mlc_snv_cache", "mlc_indel_cache", "regional_constraint_cache"
]

def load_cfg(p: str) -> dict:
    """Load and validate config.json; set defaults for optional fields."""
    with open(p, "r") as f:
        cfg = json.load(f)
    cfg.setdefault("vep_bin", "vep")
    cfg.setdefault("python", sys.executable)
    cfg.setdefault("post_script", "./hail/process_variants_combine.py")

    missing = [k for k in REQ_KEYS if not cfg.get(k)]
    if missing:
        raise SystemExit(f"ERROR: Missing keys in config.json: {missing}")

    if cfg["pipeline_mode"] not in {"population", "disease"}:
        raise SystemExit("ERROR: pipeline_mode must be 'population' or 'disease'.")

    if cfg["pipeline_mode"] == "disease" and not cfg.get("disease_meta_file"):
        raise SystemExit("ERROR: disease mode requires 'disease_meta_file'.")

    return cfg

def read_proband_list(meta_path: str) -> list[str]:
    """
    Read the disease metadata table and return a de-duplicated, sorted list
    of proband sample IDs.
    """
    meta = Path(meta_path)
    if not meta.exists():
        raise SystemExit(f"ERROR: disease_meta_file not found: {meta}")
    head = meta.read_text(errors='ignore').splitlines()[:1]
    sep = '\t' if (not head or '\t' in head[0] or ',' not in head[0]) else ','

    with meta.open(newline='', encoding='utf-8', errors='ignore') as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        cols = [c.strip() for c in (reader.fieldnames or [])]

        def find_col(cands):
            low = [c.lower() for c in cols]
            for c in cands:
                if c.lower() in low:
                    return cols[low.index(c.lower())]
            return None

        sample_col   = find_col(["Sample", "sample", "Sample_ID", "ID", "SampleID"])
        category_col = find_col(["Category", "category"])
        if not sample_col or not category_col:
            raise SystemExit(f"ERROR: Cannot find sample/category columns in {meta}. Columns: {cols}")

        probands = []
        for row in reader:
            if (row.get(category_col) or "").strip().lower() == "proband":
                s = (row.get(sample_col) or "").strip()
                if s:
                    probands.append(s)

    probands = sorted(set(probands))
    if not probands:
        raise SystemExit(f"ERROR: No samples with Category=Proband in {meta}")
    return probands

def filter_vcf_keep_if_any_proband(in_vcf: Path, out_vcf: Path, proband_names: list[str]) -> None:
    """
    Filter an *uncompressed* VCF to keep only records where at least one
    proband sample column is present.
    """
    keep = set(proband_names)
    sel_idx = None
    sel_names = None
    header_seen = False

    with Path(in_vcf).open("r", encoding="utf-8", errors="ignore") as fin, \
         Path(out_vcf).open("w", encoding="utf-8") as fout:

        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                header_seen = True
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    raise SystemExit("ERROR: VCF must contain sample columns.")
                fixed = cols[:9]
                samples = cols[9:]
                sel_names = [s for s in samples if s in keep]
                if not sel_names:
                    raise SystemExit("ERROR: None of the probands are present in the VCF header.")
                sel_idx = [9 + samples.index(s) for s in sel_names]
                fout.write("\t".join(fixed + sel_names) + "\n")
                continue

            if not header_seen:
                raise SystemExit("ERROR: Malformed VCF: data before header.")

            fields = line.rstrip("\n").split("\t")
            fixed = fields[:9]
            sample_fields = fields[9:]

            def present(val: str) -> bool:
                v = (val or "").strip()
                if v in {"", "."}:
                    return False
                first = v.split(":", 1)[0]
                if ("/" in first) or ("|" in first):
                    return first not in {"0/0", "0|0", "./."}
                return True

            keep_row = any(present(sample_fields[i - 9]) for i in sel_idx)
            if not keep_row:
                continue

            kept = [fields[i] for i in sel_idx]
            fout.write("\t".join(fixed + kept) + "\n")

    print(f"[OK] Filtered by probands. Output: {out_vcf} (samples kept: {len(sel_names)})")

# -------------------- Entry point --------------------
def main():
    ap = argparse.ArgumentParser(
        description="Run VEP and post-processing with optional disease-mode proband filtering (config-driven)."
    )
    ap.add_argument("--config", required=True, help="Path to config.json.")
    ap.add_argument("--input-vcf", required=True, help="Input *uncompressed* VCF (already gunzipped upstream).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    args = ap.parse_args()

    cfg = load_cfg(args.config)

    outdir    = Path(args.outdir).resolve()
    vep_dir   = outdir / "annotation" / "vep_vcf"
    meta_dir  = outdir / "annotation" / "metadata"
    final_dir = outdir / "annotation" / "final_outputs"
    for d in (vep_dir, meta_dir, final_dir):
        d.mkdir(parents=True, exist_ok=True)

    in_vcf = Path(args.input_vcf).resolve()
    if in_vcf.suffix == ".gz":
        raise SystemExit("ERROR: pass an *uncompressed* VCF to --input-vcf.")

    # ===== Disease mode: filter to probands first, then run VEP =====
    vep_input = in_vcf
    if cfg["pipeline_mode"] == "disease":
        probands = read_proband_list(cfg["disease_meta_file"])
        (meta_dir / "probands.list").write_text("\n".join(probands) + "\n")
        filtered_vcf = vep_dir / (in_vcf.stem + ".proband_only.filtered.vcf")
        if (not filtered_vcf.exists()) or args.overwrite:
            filter_vcf_keep_if_any_proband(in_vcf, filtered_vcf, probands)
        vep_input = filtered_vcf

    # ===== Run VEP =====
    out_vep = vep_dir / (vep_input.stem + "_vep.vcf")
    if out_vep.exists() and not args.overwrite:
        print(f"[SKIP] VEP output exists: {out_vep}")
    else:
        vep_cmd = [
            cfg["vep_bin"], "--offline","--cache",
            "--fasta", cfg["reference_genome_path"],
            "-i", str(vep_input),
            "-o", str(out_vep),
            "--dir_cache", cfg["vep_cache_dir"],
            "--no_stats",
            "--distance", "0",
            "--biotype", "--symbol", "--hgvs", "--variant_class",
            "--force_overwrite", "--vcf",
        ]
        run(vep_cmd)
    print(f"[+] VEP output: {out_vep}")

    # ===== Post-processing (Updated to pass 5 new caches) =====
    post = [
        cfg["python"], cfg["post_script"],
        "--vep-vcf", str(out_vep),
        "--final-output-dir", str(final_dir),
        "--gnomadcache",         cfg["gnomadcache"],
        "--clinvarcache",        cfg["clinvarcache"],
        "--mitomap-polycache",   cfg["mitomap_polycache"],
        "--mitomap-diseasecache",cfg["mitomap_diseasecache"],
        "--helixcache",          cfg["helixcache"],
        "--mitimpactcache",      cfg["mitimpactcache"],
        "--mitotipcache",        cfg["mitotipcache"],
        "--hmtvarcache",         cfg["hmtvarcache"],
        # NEW caches
        "--napogeecache",        cfg["napogeecache"],
        "--tapogeecache",        cfg["tapogeecache"],
        "--mlc_snv_cache",       cfg["mlc_snv_cache"],
        "--mlc_indel_cache",     cfg["mlc_indel_cache"],
        "--regional_constraint_cache", cfg["regional_constraint_cache"],
        "--pipeline-mode",       cfg["pipeline_mode"],
    ]
    if cfg["pipeline_mode"] == "disease":
        post += ["--disease-meta-file", cfg["disease_meta_file"]]

    run(post)
    print("[SUCCESS] Finished. Final outputs:", final_dir)

if __name__ == "__main__":
    main()