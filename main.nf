#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonBuilder

// ======================================================================
// CONFIGURATION & INITIALIZATION
// ======================================================================

// Print pipeline header with Work Directory to console at launch
log.info """
PIPELINE START
======================================================================
Input file:        ${params.input}
Output directory:  ${params.outdir}
Work directory:    ${workflow.workDir}
Summary log:       ${params.summary_log}
Launch directory:  ${workflow.launchDir}
======================================================================
"""

// ======================================================================
// TOP-LEVEL CHANNELS (avoid file() at operators here)
// ======================================================================
if( !params.input ) {
    System.err.println("FATAL: Missing --input parameter.")
    System.exit(1)
}
def _in = file(params.input)
if( !_in.exists() ) {
    System.err.println("FATAL: Input file does not exist: ${params.input}")
    System.exit(1)
}

ch_ref_fasta_val = Channel.value(params.wdl_inputs.ref_fasta)

ch_cromwell_conf  = Channel.value(file("${baseDir}/cromwell.conf"))
ch_hail_directory = Channel.value(file(params.hail_script))

ch_input_rows = Channel
    .fromPath(_in)
    .splitCsv(header:false, sep:'\t')
    .filter { row -> row && row.size() >= 2 && row[0] && row[0].trim() && row[1] }
    .map { row ->
        def meta = [ id: (row[0] as String) ]
        def f1 = (row[1] as String)
        def f2 = (row.size() > 2 ? (row[2] as String) : null)
                if (f1.endsWith('.cram')) {
            meta.cram = file(f1)
            meta.crai = f2 ? file(f2) : file("${f1}.crai")
        }

        def ft = f1.endsWith('.fastq.gz') || f1.endsWith('.fq.gz') ? 'FASTQ'
               : f1.endsWith('.cram') ? 'CRAM'
               : f1.endsWith('.bam')  ? 'BAM'
               : 'UNKNOWN'
        [ meta, ft, f1, f2 ]
    }
    .groupTuple()
    .map { meta, ft_list, f1_list, f2_list ->
        def unique_types = ft_list.unique()
        if( unique_types.size() > 1 )
            error "Sample ${meta.id} contains mixed input file types: ${unique_types}."
        def ft = unique_types[0]
        
        // CHECK 1: 
        def annotation_done = file("${params.outdir}/${meta.id}/annotation/final_outputs/variant_list.txt").exists()
        
        // CHECK 2: 
        def vcf = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.final.split.vcf")
        def ctm = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.haplocheck_contamination.txt")
        def cov = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.per_base_coverage.tsv")
        def wdl_done = vcf.exists() && ctm.exists() && cov.exists()

        meta.annotation_done = annotation_done
        meta.wdl_done = wdl_done

        def pairs = (0..<f1_list.size()).collect { i -> [ f1_list[i], (f2_list ? f2_list[i] : null) ] }
        [ meta, ft, pairs ]
    }
    .filter { meta, ft, pairs ->
        if (meta.annotation_done) {
            log.info "Skipping Sample: ${meta.id} (Full results already exist)"
            return false 
        }
        return true 
    }

def ch_fastqs = ch_input_rows.filter { meta, ft, pairs -> ft == 'FASTQ' }
def ch_crams  = ch_input_rows.filter { meta, ft, pairs -> ft == 'CRAM'  }
def ch_bams   = ch_input_rows.filter { meta, ft, pairs -> ft == 'BAM'   }

ch_skipped_wdl_results = ch_input_rows
    .filter { meta, ft, pairs -> meta.wdl_done }
    .map { meta, ft, pairs ->
        def input_cram = file(pairs[0][0])
        def input_crai = pairs[0][1] ? file(pairs[0][1]) : file("${pairs[0][0]}.crai")
        
        meta.cram = input_cram
        meta.crai = input_crai
        
        def vcf = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.final.split.vcf")
        def ctm = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.haplocheck_contamination.txt")
        def cov = file("${params.outdir}/${meta.id}/variant_calling/${meta.id}.per_base_coverage.tsv")
        
        return tuple(meta, vcf, ctm, cov, input_cram, input_crai)
    }

// ======================================================================
/*                              PROCESSES                               */
// ======================================================================

process ALIGN_AND_UNSORT {
  tag "BWA-MEM on ${meta.id} (Pair ${meta.pair_id})"
  label 'alignment_related'

  input:
  tuple val(meta), path(read1), path(read2)
  val ref_fasta

  output:
  tuple val(meta), path("${meta.id}_${meta.pair_id}.unsorted.bam"), emit: unsorted_bam

  script:
  def rgid = "${meta.id}_${meta.pair_id}"
  def rg   = "'@RG\\tID:${rgid}\\tSM:${meta.id}\\tPL:ILLUMINA'"
  """
  #!/usr/bin/env bash
  set -euo pipefail
  bwa mem -K 100000000 -v 3 -Y -t ${task.cpus} -R ${rg} ${ref_fasta} ${read1} ${read2} \\
    | samtools view -@ ${task.cpus} -b -o ${meta.id}_${meta.pair_id}.unsorted.bam -
  rm -f ${read1} ${read2}
  """
}

process SORT_AND_CONVERT_TO_CRAM {
  tag "Sort & CRAM for ${meta.id} (Pair ${meta.pair_id})"
  label 'alignment_related'

  input:
  tuple val(meta), path(unsorted_bam)
  val ref_fasta

  output:
  tuple val(meta), path("${meta.id}_${meta.pair_id}.cram"), path("${meta.id}_${meta.pair_id}.cram.crai"), emit: single_cram

  script:
  def ubam = unsorted_bam.getName()
  """
  set -euo pipefail
  SORTED_BAM=\"${meta.id}_${meta.pair_id}.sorted.bam\"
  OUTPUT_CRAM=\"${meta.id}_${meta.pair_id}.cram\"

  samtools sort -@ ${task.cpus} -o "\$SORTED_BAM" ${ubam}
  samtools view -@ ${task.cpus} -T ${ref_fasta} -C -o "\$OUTPUT_CRAM" "\$SORTED_BAM"
  samtools index "\$OUTPUT_CRAM"
  rm -f "\$SORTED_BAM" "${ubam}"
  """
}

process MERGE_CRAMS {
  tag "Merge CRAMs for ${meta.id}"
  label 'alignment_related'

  publishDir "${params.outdir}", mode:'copy', pattern: "*.{cram,crai}", saveAs: { fn -> "${meta.id}/alignment/${fn}" }

  input:
  tuple val(meta), path(crams), path(crais)

  output:
  tuple val(meta), path("${meta.id}.merged.cram"), path("${meta.id}.merged.cram.crai"), emit: merged_cram

  script:
  def cram_files = crams.join(' ')
  def crai_files = crais.join(' ')
  """
  #!/usr/bin/env bash
  set -euo pipefail
  OUTPUT_CRAM="${meta.id}.merged.cram"
  samtools merge -@ ${task.cpus} -f -O CRAM --output-fmt-option version=3.0 -o "\$OUTPUT_CRAM" ${cram_files}
  samtools index "\$OUTPUT_CRAM"
  rm -f ${cram_files} ${crai_files}
  """
}

process CONVERT_BAM_TO_CRAM {
  tag "Convert BAM to CRAM for ${meta.id}"
  label 'alignment_related'

  publishDir "${params.outdir}", mode:'copy', pattern:"*.{cram,crai}", saveAs: { fn -> "${meta.id}/alignment/${fn}" }

  input:
  tuple val(meta), path(input_bam)
  val ref_fasta

  output:
  tuple val(meta), path("${meta.id}.merged.cram"), path("${meta.id}.merged.cram.crai"), emit: merged_cram

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  OUTPUT_CRAM="\${meta.id}.merged.cram"
  SORTED_BAM="\${meta.id}.sorted.bam"

  samtools sort -@ ${task.cpus} -o "\$SORTED_BAM" ${input_bam}
  samtools view -@ ${task.cpus} -T ${ref_fasta} -C -o "\$OUTPUT_CRAM" "\$SORTED_BAM"
  samtools index "\$OUTPUT_CRAM"
  rm -f "\$SORTED_BAM"
  """
}

process GENERATE_CRAM_TSV {
  tag "CRAM TSV for ${meta.id}"
  label 'generation_related'

  input:
  tuple val(meta), path(cram), path(crai)

  output:
  tuple val(meta), path("${meta.id}_cram_list.tsv"), emit: tsv

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail

  # More tolerant existence check: target exists (-e) or is a symlink (-L)
  if [[ ! -e "${cram}" && ! -L "${cram}" ]]; then
    echo "[ERROR] CRAM path is neither existing nor a symlink: ${cram}" >&2
    exit 2
  fi
  if [[ ! -e "${crai}" && ! -L "${crai}" ]]; then
    echo "[ERROR] CRAI path is neither existing nor a symlink: ${crai}" >&2
    exit 3
  fi

  # Resolve absolute paths: prefer readlink -f; fall back to Python if unavailable
  abs_path() {
    local p="\$1"
    if readlink -f "\$p" >/dev/null 2>&1; then
      readlink -f "\$p"
    else
      python - <<'PY' "\$p"
import os, sys
p = sys.argv[1]
try:
    print(os.path.realpath(p))
except Exception:
    print(os.path.abspath(p))
PY
    fi
  }

  cram_path=\$(abs_path "${cram}")
  crai_path=\$(abs_path "${crai}")

  echo -e "\${cram_path}\t\${crai_path}" > ${meta.id}_cram_list.tsv
  echo "[OK] Wrote ${meta.id}_cram_list.tsv"
  """
}

process GENERATE_WDL_JSON {
  tag "WDL JSON for ${meta.id}"
  label 'generation_related'

  input:
  tuple val(meta), path(cram_tsv)

  output:
  tuple val(meta), path("${meta.id}_wdl_inputs.json"), emit: json

  script:
  // Build namespaced WDL inputs once, substituting the TSV path at runtime.
  def wdl_inputs = params.wdl_inputs.collectEntries { k,v -> ["MitochondriaMultiSamplePipeline.${k}", v] }
  wdl_inputs["MitochondriaMultiSamplePipeline.inputSamplesFile"] = "___TSV_PATH_PLACEHOLDER___"
  def template = new JsonBuilder(wdl_inputs).toPrettyString()
  """
  #!/usr/bin/env bash
  set -euo pipefail
  TSV_PATH=\$(readlink -f ${cram_tsv})
  JSON_TEMPLATE=\$(printf '%s' '${template}')
  echo "\${JSON_TEMPLATE}" | sed "s|___TSV_PATH_PLACEHOLDER___|\${TSV_PATH}|g" > ${meta.id}_wdl_inputs.json
  """
}

process RUN_WDL_VARIANT_CALLING {
  tag "Variant Calling ${meta.id}"
  label 'wdl_related'
 
  maxForks 3

  input:
  tuple val(meta), path(wdl_inputs_json)
  path cromwell_config

  output:
  tuple val(meta),
        path("nxf_emit/*.final.split.vcf"),
        path("nxf_emit/*.haplocheck_contamination.txt"),
        path("nxf_emit/*.per_base_coverage.tsv"),
        emit: wdl_files
  path("nxf_emit/*"), optional: true, emit: wdl_dump

  script: 
"""
  #!/usr/bin/env bash
  set -euo pipefail

  # --- Dynamically monitoring ---
  MAX_ALLOWED_JOBS=120
  CHECK_INTERVAL=30
  echo "[INFO] Monitoring cluster load for user \${USER}..."

  while true; do
      CURRENT_COUNT=\$(squeue -u \${USER} -h -t all | wc -l)
      
      if (( CURRENT_COUNT < MAX_ALLOWED_JOBS )); then
          echo "[OK] Cluster load is \${CURRENT_COUNT}. Proceeding with ${meta.id}..."
          break
      else
          echo "[WAIT] High load (\${CURRENT_COUNT} jobs). Limit is \${MAX_ALLOWED_JOBS}. Sleeping \${CHECK_INTERVAL}s..."
          sleep \${CHECK_INTERVAL}
      fi
  done
  # ---------------------------------------------

  TARGET_DIR="${params.outdir}/${meta.id}/variant_calling"
  mkdir -p "\${TARGET_DIR}"

  cat > cromwell_options.json <<EOF
  {
    "final_workflow_outputs_dir": "\${TARGET_DIR}",
    "default_runtime_attributes": {
        "queue": "${params.cromwell_options.queue ?: ''}",
        "cpus": ${params.cromwell_options.cpus},
        "memory": ${params.cromwell_options.memory},
        "runtime_minutes": ${params.cromwell_options.runtime_minutes}
      }
  }
EOF
 cat > sbatch_throttle.sh <<'EOS'
    #!/usr/bin/env bash
    set -euo pipefail
    PER_HOUR=50
    MIN_GAP=\$(( 3600 / PER_HOUR ))
    STATE_DIR="\${HOME}/.sbatch_rate"
    LOCK_FILE="\${STATE_DIR}/lock"
    TS_FILE="\${STATE_DIR}/last_submit.ts"
    mkdir -p "\${STATE_DIR}"
    exec 200>"\${LOCK_FILE}"
    flock 200
    now=\$(date +%s); last=0
    [[ -f "\${TS_FILE}" ]] && read -r last < "\${TS_FILE}" || true
    delta=\$(( now - last ))
    if (( delta < MIN_GAP )); then
      sleep \$(( MIN_GAP - delta ))
    fi
    date +%s > "\${TS_FILE}"
    unset SLURM_CONF || true
    exec sbatch "\$@"
EOS
    chmod +x sbatch_throttle.sh

  export PATH="\$PWD:\$PATH"

  java -Dconfig.file=${cromwell_config} -jar ${params.cromwell_jar} run ${params.wdl_script} --inputs ${wdl_inputs_json} --options cromwell_options.json

  # Flatten important outputs from deep cromwell dirs into TARGET_DIR
  shopt -s nullglob
  patterns=(
    "*.final.split.vcf" "*.final.split.vcf.idx"
    "*.final.vcf.gz" "*.final.vcf.gz.tbi"
    "*.vcf.gz" "*.vcf.gz.tbi"
    "*.vcf.gz.numt.vcf.gz" "*.vcf.gz.numt.vcf.gz.tbi"
    "*.splitAndPassOnly.vcf" "*.splitAndPassOnly.vcf.idx"
    "*.per_base_coverage.tsv"
    "*.haplocheck_contamination.txt"
    "*.metrics" "metrics.txt" "theoretical_sensitivity.txt"

  )

  for pat in "\${patterns[@]}"; do
    while IFS= read -r -d '' f; do
      base="\${f##*/}"
      if [[ ! -e "\${TARGET_DIR}/\${base}" ]] || ! cmp -s "\$f" "\${TARGET_DIR}/\${base}"; then
        cp -fL "\$f" "\${TARGET_DIR}/\${base}"
      fi
    done < <(find -L "\${TARGET_DIR}" -mindepth 2 -type f -name "\$pat" -print0)
  done

  WORK_EXEC="\$PWD/cromwell-executions"
  if [[ -d "\${WORK_EXEC}" ]]; then
    for pat in "\${patterns[@]}"; do
      while IFS= read -r -d '' f; do
        base="\${f##*/}"
        if [[ ! -e "\${TARGET_DIR}/\${base}" ]] || ! cmp -s "\$f" "\${TARGET_DIR}/\${base}"; then
          cp -fL "\$f" "\${TARGET_DIR}/\${base}"
        fi
      done < <(find -L "\${WORK_EXEC}" -mindepth 2 -type f -name "\$pat" -print0)
    done
  fi

  # Preserve regular-mt realigned outputs with canonical sample-level names
  copy_regular_realigned() {
    local search_root="\$1"
    local regular_bam="${meta.id}.realigned.bam"
    local regular_bai="${meta.id}.realigned.bai"
    [[ -d "\${search_root}" ]] || return 0

    local bam_src=""
    bam_src=\$(find -L "\${search_root}" -type f -path '*/call-AlignToMt/*' -name '*.realigned.bam' -print | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bam_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${regular_bam}" ]] || ! cmp -s "\${bam_src}" "\${TARGET_DIR}/\${regular_bam}"; then
        cp -fL "\${bam_src}" "\${TARGET_DIR}/\${regular_bam}"
      fi
    fi

    local bai_src=""
    bai_src=\$( ( find -L "\${search_root}" -type f -path '*/call-AlignToMt/*' -name '*.realigned.bai' -print; find -L "\${search_root}" -type f -path '*/call-AlignToMt/*' -name '*.realigned.bam.bai' -print ) | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bai_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${regular_bai}" ]] || ! cmp -s "\${bai_src}" "\${TARGET_DIR}/\${regular_bai}"; then
        cp -fL "\${bai_src}" "\${TARGET_DIR}/\${regular_bai}"
      fi
    fi
  }

  # Preserve shifted-mt realigned outputs with canonical sample-level names
  copy_shifted_realigned() {
    local search_root="\$1"
    local shifted_bam="${meta.id}.realigned.shifted.bam"
    local shifted_bai="${meta.id}.realigned.shifted.bai"
    [[ -d "\${search_root}" ]] || return 0

    local bam_src=""
    bam_src=\$(find -L "\${search_root}" -type f -path '*/call-AlignToShiftedMt/*' -name '*.realigned.bam' -print | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bam_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${shifted_bam}" ]] || ! cmp -s "\${bam_src}" "\${TARGET_DIR}/\${shifted_bam}"; then
        cp -fL "\${bam_src}" "\${TARGET_DIR}/\${shifted_bam}"
      fi
    fi

    local bai_src=""
    bai_src=\$( ( find -L "\${search_root}" -type f -path '*/call-AlignToShiftedMt/*' -name '*.realigned.bai' -print; find -L "\${search_root}" -type f -path '*/call-AlignToShiftedMt/*' -name '*.realigned.bam.bai' -print ) | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bai_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${shifted_bai}" ]] || ! cmp -s "\${bai_src}" "\${TARGET_DIR}/\${shifted_bai}"; then
        cp -fL "\${bai_src}" "\${TARGET_DIR}/\${shifted_bai}"
      fi
    fi
  }

  copy_regular_realigned "\${TARGET_DIR}"
  copy_regular_realigned "\${WORK_EXEC}"
  copy_shifted_realigned "\${TARGET_DIR}"
  copy_shifted_realigned "\${WORK_EXEC}"

  # Keep canonical sample-level BAM/BAI names (e.g. sample.bam / sample.bai)
  # from SubsetBamToChrM outputs.
  copy_sample_level_bam() {
    local search_root="\$1"
    local sample_bam="${meta.id}.bam"
    local sample_bai="${meta.id}.bai"
    [[ -d "\${search_root}" ]] || return 0

    local bam_src=""
    bam_src=\$(find -L "\${search_root}" -type f -path '*/call-SubsetBamToChrM/*' -name '*.bam' -print | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bam_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${sample_bam}" ]] || ! cmp -s "\${bam_src}" "\${TARGET_DIR}/\${sample_bam}"; then
        cp -fL "\${bam_src}" "\${TARGET_DIR}/\${sample_bam}"
      fi
    fi

    local bai_src=""
    bai_src=\$( ( find -L "\${search_root}" -type f -path '*/call-SubsetBamToChrM/*' -name '*.bai' -print; find -L "\${search_root}" -type f -path '*/call-SubsetBamToChrM/*' -name '*.bam.bai' -print ) | LC_ALL=C sort | head -n1 || true)
    if [[ -n "\${bai_src}" ]]; then
      if [[ ! -e "\${TARGET_DIR}/\${sample_bai}" ]] || ! cmp -s "\${bai_src}" "\${TARGET_DIR}/\${sample_bai}"; then
        cp -fL "\${bai_src}" "\${TARGET_DIR}/\${sample_bai}"
      fi
    fi
  }

  copy_sample_level_bam "\${TARGET_DIR}"
  copy_sample_level_bam "\${WORK_EXEC}"

  # Step 3: remove undesired BAM/BAI outputs and keep only canonical set
  keep_bam_names=(
    "${meta.id}.bam"
    "${meta.id}.bai"
    "${meta.id}.realigned.bam"
    "${meta.id}.realigned.bai"
    "${meta.id}.realigned.shifted.bam"
    "${meta.id}.realigned.shifted.bai"
  )

  for f in "\${TARGET_DIR}"/*.bam "\${TARGET_DIR}"/*.bai "\${TARGET_DIR}"/*.bam.bai; do
    [[ -e "\$f" ]] || continue
    base="\${f##*/}"
    keep=false
    for k in "\${keep_bam_names[@]}"; do
      if [[ "\${base}" == "\${k}" ]]; then
        keep=true
        break
      fi
    done
    if [[ "\${keep}" == false ]]; then
      rm -f "\$f"
    fi
  done

  # Existence checks for three core artifacts
  VCF=\$(ls -1 "\${TARGET_DIR}"/*.final.split.vcf 2>/dev/null | head -n1 || true)
  CTM=\$(ls -1 "\${TARGET_DIR}"/*.haplocheck_contamination.txt 2>/dev/null | head -n1 || true)
  COV=\$(ls -1 "\${TARGET_DIR}"/*.per_base_coverage.tsv 2>/dev/null | head -n1 || true)

  if [[ -z "\${VCF}" || -z "\${CTM}" || -z "\${COV}" ]]; then
    echo "[ERROR] Missing required WDL outputs for ${meta.id}:" >&2
    echo "  VCF: \${VCF:-<none>}"  >&2
    echo "  CTM: \${CTM:-<none>}"  >&2
    echo "  COV: \${COV:-<none>}"  >&2
    exit 1
  fi

  # Emit symlinks for Nextflow outputs
  EMIT_DIR="nxf_emit"
  rm -rf "\${EMIT_DIR}"; mkdir -p "\${EMIT_DIR}"
  ln -s "\${VCF}" "\${EMIT_DIR}/\${VCF##*/}"
  ln -s "\${CTM}" "\${EMIT_DIR}/\${CTM##*/}"
  ln -s "\${COV}" "\${EMIT_DIR}/\${COV##*/}"

  for pat in "\${patterns[@]}"; do
    for f in "\${TARGET_DIR}"/\$pat; do
      [[ -e "\$f" ]] || continue
      base="\${f##*/}"
      [[ -e "\${EMIT_DIR}/\${base}" ]] || ln -s "\$f" "\${EMIT_DIR}/\${base}"
    done
  done

  echo "[OK] Final variant_calling outputs for ${meta.id}:"
  ls -lh "\${TARGET_DIR}" || true

  echo "[OK] nxf_emit outputs for ${meta.id}:"
  ls -lh "\${EMIT_DIR}" || true
"""
}

process CALCULATE_MTCN {
  tag "Calculate_mtCN for ${meta.id}"
  label 'mtCN_related'

  publishDir "${params.outdir}/${meta.id}/mtCN", mode: 'copy'

  input:
  tuple val(meta), path(cram), path(crai), path(mt_coverage), path(haplocheck_file)
  path ref_fasta
  path hail_dir

  output:
  tuple val(meta), path("mtCN_summary.txt"), emit: mtcn_summary

  script:
  def intervals_arg = params.wgs_intervals ? "--intervals ${params.wgs_intervals}" : ""
  """
  #!/usr/bin/env bash
  set -euo pipefail
  mkdir -p mtcn_out

  REF_ABS=\$(readlink -f ${ref_fasta})
  CRAM_ABS=\$(readlink -f ${cram})
  CRAI_ABS=\$(readlink -f ${crai})
  MT_COV_ABS=\$(readlink -f ${mt_coverage})

  python ${hail_dir}/calculate_mtcn_mosdepth.py \\
      --cram "\${CRAM_ABS}" \\
      --crai "\${CRAI_ABS}" \\
      --ref_fasta "\${REF_ABS}" \\
      --mt_coverage "\${MT_COV_ABS}" \\
      --haplocheck ${haplocheck_file} \\
      --output mtcn_out \\
      --mosdepth ${params.mosdepth} \\
      --min_mt_cov ${params.min_mt_cov ?: 100} \\
      --min_mtcn ${params.min_mtcn ?: 50} \\
      --max_mtcn ${params.max_mtcn ?: 500} \\
      --max_contam ${params.max_contam ?: 0.02}

  mv mtcn_out/mtCN_summary.txt ./mtCN_summary.txt
  echo "[OK] mtCN summary generated for ${meta.id}"
  """
}

process SAMPLE_LEVEL_FILTER {
  tag "Sample filter ${meta.id}"
  label 'mtCN_related'

  publishDir "${params.outdir}/${meta.id}", mode: 'copy'

  input:
  tuple val(meta), path(mtcn_summary)

  output:
  tuple val(meta), path(".pass.ok"), optional: true, emit: pass_signal
  path "${meta.id}.individual_entry.txt", emit: log_entry

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail

  results=\$(python - "${mtcn_summary}" <<'PY'
import sys

try:
    summary_path = sys.argv[1]
    with open(summary_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        if len(lines) < 2:
            raise ValueError("Summary file is empty or missing data row.")
            
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        res_dict = dict(zip(header, data))
        
        mt_cov = res_dict.get('Mean_mtDNA_Coverage', 'N/A')
        mtcn   = res_dict.get('mtCN_final', 'N/A')
        contam_status = res_dict.get('Contamination', 'NO').upper()
        
        orig_pass = res_dict.get('Pass_Filter', 'False').lower() in ['true', 'pass', '1']
        
        final_pass = "TRUE" if (orig_pass and contam_status != "YES") else "FALSE"
        
        print(f"{final_pass}|{mt_cov}|{mtcn}|{contam_status}")
except Exception as e:
    sys.stderr.write(f"Python Error: {str(e)}\\n")
    print("ERROR|N/A|N/A|N/A")
PY
)

  IFS='|' read -r PASS_FLAG COV MTCN CONTAM <<< "\${results}"

  if [[ "\${PASS_FLAG}" == "TRUE" ]]; then
    echo -e "${meta.id}\tSUCCESS\tTRUE\t\${COV}\t\${MTCN}\t\${CONTAM}\tQC Pass" > "${meta.id}.individual_entry.txt"
    touch .pass.ok
  else
    echo -e "${meta.id}\tSUCCESS\tFALSE\t\${COV}\t\${MTCN}\t\${CONTAM}\tQC Failed or Contaminated" > "${meta.id}.individual_entry.txt"
  fi
  """
}

process ANNOTATE_INDIVIDUAL_VCF {
  tag "Annotate VCF ${meta.id}"
  label 'vep_related'

  publishDir "${params.outdir}", mode:'copy', saveAs: { fn -> "${meta.id}/${fn}" }

  input:
  tuple val(meta), path(vcf), path(contamination), path(coverage)
  path hail_dir

  output:
  tuple val(meta), path("annotation/final_outputs/variant_list_prefiltering.txt"), emit: prefilter_out

  script:
  def cfg = new JsonBuilder( params.hail_pipeline_config ).toPrettyString()
  def in_dir = "wdl_outputs"
  """
  #!/usr/bin/env bash
  set -euo pipefail

  printf '%s' '${cfg}' > config.json
  rm -rf annotation
  mkdir -p ${in_dir}/vcfs ${in_dir}/contamination ${in_dir}/coverage
  mkdir -p annotation/{vep_vcf,metadata,final_outputs}
  mkdir -p ${params.outdir}/merged_results/wdl_output

  cp "${vcf}"           "${in_dir}/vcfs/${meta.id}.merged.final.split.vcf"
  cp "${contamination}" "${in_dir}/contamination/${meta.id}.merged.haplocheck_contamination.txt"
  cp "${coverage}"      "${in_dir}/coverage/${meta.id}.merged.per_base_coverage.tsv"

  cp "${vcf}"           "${params.outdir}/merged_results/wdl_output/${meta.id}.merged.final.split.vcf"
  cp "${contamination}" "${params.outdir}/merged_results/wdl_output/${meta.id}.merged.haplocheck_contamination.txt"
  cp "${coverage}"      "${params.outdir}/merged_results/wdl_output/${meta.id}.merged.per_base_coverage.tsv"

  python ${hail_dir}/add_annotation_single.py --config config.json

  FINAL_OUTPUT_FILE="annotation/final_outputs/variant_list_prefiltering.txt" 
  
  if [[ -s "\${FINAL_OUTPUT_FILE}" ]]; then
    echo "[SUCCESS] Annotation OK for ${meta.id} -> \${FINAL_OUTPUT_FILE}"
    touch .annotate_complete
  else
    echo "ERROR: Expected final output '\${FINAL_OUTPUT_FILE}' not found or empty." >&2
    exit 1
  fi
  """
}

// ======================================================================
// WORKFLOW
// ======================================================================
workflow {

  // ---------------- FASTQ path ----------------
  ch_fastq_pairs = ch_fastqs.flatMap { meta, type, pairs ->
    pairs.indexed().collect { i, p -> [ meta + [pair_id: i], p[0], p[1] ] }
  }

  ALIGN_AND_UNSORT(ch_fastq_pairs, ch_ref_fasta_val)
  SORT_AND_CONVERT_TO_CRAM(ALIGN_AND_UNSORT.out.unsorted_bam, ch_ref_fasta_val)

  ch_crams_from_fastq =
    SORT_AND_CONVERT_TO_CRAM.out.single_cram
      .map { meta, cram, crai -> tuple(meta.id, tuple(meta, cram, crai)) }
      .groupTuple()
      .map { sid, tuples ->
        def meta  = tuples.first()[0]
        def crams = tuples.collect { it[1] }
        def crais = tuples.collect { it[2] }
        tuple(meta, crams, crais)
      }

  ch_fastq_split = ch_crams_from_fastq.branch {
    merge : it[1].size() > 1
    single: it[1].size() <= 1
  }
  ch_merge_inputs_fastq = ch_fastq_split.merge
  ch_single_fastq = ch_fastq_split.single.map { rec ->
    def meta  = rec[0]; def crams = rec[1]; def crais = rec[2]
    tuple(meta, crams[0], crais[0])
  }

  // ---------------- BAM path ----------------
  ch_bam_crams = ch_bams.flatMap { meta, type, pairs ->
    pairs.indexed().collect { i, p -> [meta + [pair_id: i], p[0]] }
  }
  CONVERT_BAM_TO_CRAM(ch_bam_crams, ch_ref_fasta_val)

  ch_final_from_bam = CONVERT_BAM_TO_CRAM.out.merged_cram
    .map { meta, cram, crai -> tuple(meta.id, tuple(meta, cram, crai)) }
    .groupTuple()
    .map { sid, tuples ->
      def meta  = tuples.first()[0]
      def crams = tuples.collect { it[1] }
      def crais = tuples.collect { it[2] }
      tuple(meta, crams, crais)
    }
    .branch {
      merge : it[1].size() > 1
      single: it[1].size() <= 1
    }

  ch_merge_inputs_bam = ch_final_from_bam.merge
  ch_single_bam = ch_final_from_bam.single.map { rec ->
    def meta  = rec[0]; def crams = rec[1]; def crais = rec[2]
    tuple(meta, crams[0], crais[0])
  }

  // ---------------- CRAM path ----------------
  ch_final_from_cram = ch_crams
    .map { meta, type, pairs -> [ meta, pairs.collect { it[0] }, pairs.collect { it[1] } ] }
    .branch {
      merge : it[1].size() > 1
      single: it[1].size() <= 1
    }

  ch_merge_inputs_cram = ch_final_from_cram.merge
  ch_single_cram = ch_final_from_cram.single.map { rec ->
    def meta  = rec[0]; def crams = rec[1]; def crais = rec[2]
    tuple(meta, crams[0], crais[0])
  }

  // Merge-only branch across sources (FASTQ/BAM/CRAM)
  ch_all_merge_inputs = ch_merge_inputs_fastq.mix(ch_merge_inputs_bam).mix(ch_merge_inputs_cram)
  MERGE_CRAMS(ch_all_merge_inputs)
  ch_merged_results = MERGE_CRAMS.out.merged_cram  // (meta, merged.cram, merged.crai)

  // Single-file (no merge needed) branch
  ch_all_singles = ch_single_fastq.mix(ch_single_bam).mix(ch_single_cram)

  // Unified (meta, cram, crai) for downstream
 ch_all_crams = ch_merged_results.mix(ch_all_singles)

    // 1. WDL 门控逻辑
    ch_wdl_gate = ch_all_crams
        .map { meta, cram, crai ->
            def out_base = new File(params.outdir).absolutePath
            def dir = new File("${out_base}/${meta.id}/variant_calling")
            def files = (dir.exists() ? dir.listFiles() : null)
            def vcfP = files?.find { it.isFile() && it.name.endsWith(".final.split.vcf") }?.toString()
            def ctmP = files?.find { it.isFile() && it.name.endsWith(".haplocheck_contamination.txt") }?.toString()
            def covP = files?.find { it.isFile() && it.name.endsWith(".per_base_coverage.tsv") }?.toString()
            def wdl_exists = (vcfP && ctmP && covP)
            
            def updated_meta = meta + [wdl_exists: wdl_exists, vcf_path: vcfP, ctm_path: ctmP, cov_path: covP]
            return [updated_meta, cram, crai]
        }
        .branch {
            to_run: !it[0].wdl_exists
            exists:  it[0].wdl_exists
        }

    // 准备 WDL 运行
    GENERATE_CRAM_TSV(ch_wdl_gate.to_run)
    GENERATE_WDL_JSON(GENERATE_CRAM_TSV.out.tsv)
    RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON.out.json, ch_cromwell_conf)

    // 汇合 WDL 结果 (使用 .mix 产生 Queue Channel)
    ch_existing_wdl = ch_wdl_gate.exists.map { m, cr, ci -> 
        [ m, file(m.vcf_path), file(m.ctm_path), file(m.cov_path) ] 
    }
    
    // 重要：使用 .broadcast() 解决多下游消费问题
  ch_wdl_combined = RUN_WDL_VARIANT_CALLING.out.wdl_files
        .mix(ch_existing_wdl)
        .multiMap { meta, vcf, ctm, cov ->
            for_mtcn: [ meta, vcf, ctm, cov ]
            for_ann:  [ meta, vcf, ctm, cov ]
        }

    // 2. mtCN 准备 (消费 ch_wdl_split.for_mtcn)
    // 2. mtCN 准备
    ch_mtcn_prep = ch_wdl_combined.for_mtcn.map { meta, vcf, ctm, cov ->
        def out_base = new File(params.outdir).absolutePath
        def mtcn_summary = new File("${out_base}/${meta.id}/mtCN/mtCN_summary.txt")
        def has_mtcn = mtcn_summary.exists() && mtcn_summary.length() > 0
        
        // 健壮性检查：如果 output 目录没文件，尝试从 meta 中取原始路径
        def out_cram_file = file("${params.outdir}/${meta.id}/alignment/${meta.id}.merged.cram")
        def out_crai_file = file("${params.outdir}/${meta.id}/alignment/${meta.id}.merged.crai")
        
        def final_cram = out_cram_file.exists() ? out_cram_file : file(meta.cram)
        def final_crai = out_crai_file.exists() ? out_crai_file : file(meta.crai)
        
        // 返回 5 个元素的元组：meta, cram, crai, cov, ctm
        return [ meta + [mtcn_done: has_mtcn], final_cram, final_crai, cov, ctm ]
    }

    // 修复点：显式为 ch_mtcn_split 赋值
    ch_mtcn_split = ch_mtcn_prep.branch {
        to_run: !it[0].mtcn_done
        exists: it[0].mtcn_done
    }

    // 运行计算
    def ch_to_calculate = ch_mtcn_split.to_run
    CALCULATE_MTCN(ch_to_calculate, ch_ref_fasta_val, ch_hail_directory)

    // 3. 统一 mtCN 结果输出 (确保 downstream 拿到的是计算产出的 Path 对象)
    ch_existing_mtcn = ch_mtcn_split.exists.map { m, cr, ci, cov, ctm -> 
    [ m, file("${params.outdir}/${m.id}/mtCN/mtCN_summary.txt") ] 
}
    
    // 关键：这里直接消费 CALCULATE_MTCN.out，保证了 Fresh Run 时的线性依赖
    ch_all_mtcn_results = CALCULATE_MTCN.out.mtcn_summary
        .mix(ch_existing_mtcn)

    // 4. 运行 QC 过滤
    SAMPLE_LEVEL_FILTER(ch_all_mtcn_results)

    // 5. 最终汇合 (使用 ID 作为 Key 进行 join)
    // 5. 最终汇合
    ch_qc_pass = SAMPLE_LEVEL_FILTER.out.pass_signal
        .map { meta, ok -> [ meta.id.toString(), meta ] }

    // 使用 multiMap 的 for_ann 分支
    ch_wdl_for_ann = ch_wdl_combined.for_ann
        .map { meta, vcf, ctm, cov -> [ meta.id.toString(), vcf, ctm, cov ] }

    ch_annotate_inputs = ch_qc_pass
        .join(ch_wdl_for_ann)
        .map { sid, meta, vcf, ctm, cov ->
            log.info "[PROCESS] Starting Annotation for Sample: ${sid}"
            return [ meta, vcf, ctm, cov ]
        }

    // 6. 运行最终步骤
    ANNOTATE_INDIVIDUAL_VCF(ch_annotate_inputs, ch_hail_directory)
}

// ======================================================================
// WORKFLOW HOOKS
// ======================================================================
workflow.onComplete {
    if( workflow.success ) {
        log.info """
        Pipeline completed successfully. 
        Output files : ${params.outdir}
        """.stripIndent()
    }
}

workflow.onError { e ->
  log.error """
----------------------------------------------------
ERROR: Pipeline execution failed!
Message: ${e?.message}
----------------------------------------------------
"""
}
