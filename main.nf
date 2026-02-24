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

def ch_fastqs = ch_input_rows.filter { meta, ft, pairs -> ft == 'FASTQ' && !meta.wdl_done }
def ch_crams  = ch_input_rows.filter { meta, ft, pairs -> ft == 'CRAM'  && !meta.wdl_done }
def ch_bams   = ch_input_rows.filter { meta, ft, pairs -> ft == 'BAM'   && !meta.wdl_done }

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
    PER_HOUR="\${SUBMITS_PER_HOUR:-90}"
    (( PER_HOUR > 0 )) || PER_HOUR=90
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

  export SUBMITS_PER_HOUR="${params.cromwell_submit_rate_limit ?: '90'}"

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
    "*.realigned.bam" "*.realigned.bai" "*.realigned.bwa.stderr.log"
    "*.metrics" "metrics.txt" "theoretical_sensitivity.txt"
    "*.bai"
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
      base="\${f##*/}"
      [[ -e "\${EMIT_DIR}/\${base}" ]] || ln -s "\$f" "\${EMIT_DIR}/\${base}"
    done
  done
"""
}

process CALCULATE_MTCN {
  tag "Calculate_mtCN for ${meta.id}"
  label 'mtCN_related'

  publishDir "${params.outdir}/${meta.id}/mtCN", mode: 'copy'

  input:
  tuple val(meta), path(cram), path(crai), path(mt_coverage)
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

  # 使用原生 Python 提高兼容性，提取数值和 Pass 状态
  results=\$(python - "${mtcn_summary}" <<'PY'
import sys, re

try:
    with open(sys.argv[1], 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        
        # 将列名存入字典方便查找
        res_dict = dict(zip(header, data))
        
        # 匹配 calculate_mtcn.py 的原始列名
        mt_cov = res_dict.get('Mean_mtDNA_Coverage', 'N/A')
        mtcn = res_dict.get('mtCN_final', 'N/A')
        contam = res_dict.get('Contamination', 'N/A')
        # 兼容 True/TRUE/pass 等写法
        pass_val = "TRUE" if res_dict.get('Pass_Filter', 'False').lower() in ['true', 'pass', '1'] else "FALSE"
        
        print(f"{pass_val}|{mt_cov}|{mtcn}|{contam}")
except Exception as e:
    sys.stderr.write(str(e))
    print("ERROR|N/A|N/A|N/A")
PY
)

  # 解析 Python 输出
  IFS='|' read -r PASS_FLAG COV MTCN CONTAM <<< "\${results}"

  # 写入临时日志条目 (修正2：确保文件名与 output 部分定义的一致)
  if [[ "\${PASS_FLAG}" == "TRUE" ]]; then
    echo -e "${meta.id}\tSUCCESS\t\${PASS_FLAG}\t\${COV}\t\${MTCN}\t\${CONTAM}\tQC Done" > "${meta.id}.individual_entry.txt"
    touch .pass.ok
  else
    echo -e "${meta.id}\tSUCCESS\tFALSE\t\${COV}\t\${MTCN}\t\${CONTAM}\tQC Failed" > "${meta.id}.individual_entry.txt"
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
  tuple val(meta), path("annotation/**"), path(".annotate_complete"), emit: annotated_results

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

  FINAL_OUTPUT_FILE="annotation/final_outputs/variant_list.txt" 
  if [[ -s "\${FINAL_OUTPUT_FILE}" ]]; then
    echo "[SUCCESS] Annotation OK for ${meta.id} -> \${FINAL_OUTPUT_FILE}"
    touch .annotate_complete
  else
    echo "ERROR: Expected final output '\${FINAL_OUTPUT_FILE}' not found or empty." >&2
    ls -R annotation >&2 || true
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
    .map { meta, cram, crai ->
        meta.cram = cram
        meta.crai = crai
        return [ meta, cram, crai ]
    }
  // ---------------- downstream ----------------
  GENERATE_CRAM_TSV(ch_all_crams)
  GENERATE_WDL_JSON(GENERATE_CRAM_TSV.out.tsv)
  RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON.out.json, ch_cromwell_conf)

  ch_new_wdl_results = RUN_WDL_VARIANT_CALLING.out.wdl_files
    .map { meta, vcf, ctm, cov ->
        def out_cram = file("${params.outdir}/${meta.id}/alignment/${meta.id}.merged.cram")
        def out_crai = file("${params.outdir}/${meta.id}/alignment/${meta.id}.merged.cram.crai")

        def final_cram = out_cram.exists() ? out_cram : (meta.cram ?: null)
        def final_crai = out_crai.exists() ? out_crai : (meta.crai ?: null)

        if (final_cram == null) {
            error "Sample ${meta.id} has no Cram file. Check if 'cram' is defined in input TSV or alignment was skipped."
        }

        return tuple(meta, vcf, ctm, cov, final_cram, final_crai)
    }

  ch_all_wdl_final = ch_new_wdl_results.mix(ch_skipped_wdl_results)

  ch_mtcn_inputs = ch_all_wdl_final.map { meta, vcf, ctm, cov, cram, crai ->
      tuple(meta, cram, crai, cov)
  }
  
  CALCULATE_MTCN(ch_mtcn_inputs, ch_ref_fasta_val, ch_hail_directory)
  SAMPLE_LEVEL_FILTER(CALCULATE_MTCN.out.mtcn_summary)

  ch_wdl_all_keyed = ch_all_wdl_final.map { meta, vcf, ctm, cov, cram, crai ->
      tuple(meta.id as String, [meta, vcf, ctm, cov]) 
  }

  ch_annotate_inputs = SAMPLE_LEVEL_FILTER.out.pass_signal
    .map { meta, ok_file -> tuple(meta.id as String, meta) } 
    .join(ch_wdl_all_keyed)                                     
    .map { sid, meta, wdl_data -> 
        tuple(meta, wdl_data[1], wdl_data[2], wdl_data[3]) 
    }

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