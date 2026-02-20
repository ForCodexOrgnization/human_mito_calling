#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ======================================================================
// FINALIZER: merge per-sample results and produce combined outputs
// ======================================================================
log.info """
STAGE 2: FINALIZER (Data merging only)
============================================================
Merged inputs: ${params.merged_dir}
Hail scripts:  ${params.hail_script}
============================================================
"""

// ======================================================================
// INPUT DISCOVERY
// ======================================================================

def dir = file(params.merged_dir)
if ( !dir.exists() )
    error "Merged directory does not exist: ${params.merged_dir}"

def all_files = dir.listFiles()

// Select expected inputs
def vcf_files = all_files.findAll { it.name ==~ /.*merged\.final\.split\.vcf(\.gz)?$/ }
def cov_files = all_files.findAll { it.name ==~ /.*merged\.per_base_coverage\.tsv(\.gz)?$/ }

if ( vcf_files.isEmpty() )
    error "No merged.final.split.vcf(.gz) files found under ${params.merged_dir}"

if ( cov_files.isEmpty() )
    error "No merged.per_base_coverage.tsv(.gz) files found under ${params.merged_dir}"

ch_summary_entries = Channel.fromPath("${params.outdir}/*/*.individual_entry.txt").collect()

// Materialize channels
ch_vcfs     = Channel.value(vcf_files)
ch_cov_tsvs = Channel.value(cov_files)
ch_hail_dir = Channel.fromPath(params.hail_script)


// ======================================================================
// PROCESS 1: Convert per-base coverage TSV -> Hail Table/Matrix
// ======================================================================
process ANNOTATE_COVERAGE {
    tag "AnnotateCoverage_MT"
    label 'hail_related'
    publishDir "${params.outdir}/merged_results/coverage_mt", mode: 'copy'

    input:
    path coverage_tsvs
    path hail_dir

    output:
    path "coverage.mt", emit: coverage_mt
    path "coverage.ht", optional: true
    path "coverage_files.list", emit: coverage_list

    script:
    """
    set -euo pipefail

    > coverage_files.list
    for filename in ${coverage_tsvs}; do
      [[ -f "\$filename" ]] || continue
      sample_id=\$(basename "\$filename" | sed -E 's/(\\.final\\.split)?\\.per_base_coverage\\.tsv//')
      echo -e "\${sample_id}\t\$(realpath "\$filename")" >> coverage_files.list
    done

    if [[ ! -s coverage_files.list ]]; then
      echo "ERROR: coverage_files.list is empty" >&2
      exit 2
    fi

    ${params.hail_python} ${hail_dir}/annotate_coverage.py \\
      -i coverage_files.list \\
      -o coverage.ht \\
      --overwrite
    """
}


// ======================================================================
// PROCESS 2: Combine VCFs + coverage MT -> combined MT/VCF
// ======================================================================
process COMBINE_VCFS {
  tag "COMBINE_VCFS"
  label 'hail_related'
  publishDir "${params.outdir}/merged_results/combined_mt", mode: 'copy'

  input:
  path vcfs
  path coverage_mt
  path hail_dir

  output:
  path "combined_final.mt",     emit: combined_mt
  path "combined_final.vcf.bgz", emit: combined_vcf
  path "combined_final.vcf.bgz.tbi", optional: true

  script:
  """
  set -euo pipefail

  > vcf_files.list
  for filename in ${vcfs}; do
    [[ -f "\$filename" ]] || continue
    sample_id=\$(basename "\$filename" | sed -E 's/\\.merged\\.final\\.split\\.vcf(\\.gz)?\$//')
    abs_path=\$(realpath "\$filename" || echo "\$PWD/\$filename")
    echo -e "\${sample_id}\t\${abs_path}" >> vcf_files.list
  done

  if [[ ! -s vcf_files.list ]]; then
    echo "[ERROR] No VCF files collected; check file patterns." >&2
    exit 2
  fi

  ${params.hail_python} ${hail_dir}/combine_vcfs.py \\
    -i vcf_files.list \\
    -c ${coverage_mt} \\
    -a ${params.wdl_inputs.blacklisted_sites} \\
    -o . \\
    --file_suffix final \\
    --overwrite
  """
}


// ======================================================================
// PROCESS 3: Decompress combined VCF and run downstream annotation
// ======================================================================
process REFINE_ANNOTATION {
  tag "REFINE_ANNOTATION"
  label 'vep_related'
  publishDir "${params.outdir}/merged_results", mode: 'copy'

  errorStrategy 'ignore'

  input:
  path combined_vcf
  path hail_dir

  output:
  path "annotation", emit: refined_results

  script:
  // Serialize Hail/annotation configuration for the Python driver
  def cfg = params.hail_pipeline_config ?: [:]
  cfg.pipeline_mode = (params.pipeline_mode ?: 'population')  // 'population' | 'disease'
  if (cfg.pipeline_mode == 'disease' && !cfg.disease_meta_file)
      cfg.disease_meta_file = params.disease_meta_file ?: ''
  def cfgJson = new groovy.json.JsonBuilder(cfg).toPrettyString()

  """
  set -euo pipefail

  gunzip -c ${combined_vcf} > combined_final.vcf

  cat > config.json <<'JSON'
${cfgJson}
JSON

  rm -rf annotation
  mkdir -p annotation

  python ${hail_dir}/add_annotation_combine.py \\
    --config config.json \\
    --input-vcf combined_final.vcf \\
    --outdir ./ \\
    --overwrite

  """
}

process COMBINE_SUMMARY_LOGS {
    label 'hail_related'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path entry_files
    val status
    val expected_count

    output:
    path "pipeline_summary_log.txt"

    script:
    def header = """# ====================================================
# PIPELINE RUN SUMMARY
# Start Time:   ${workflow.start}
# Params:       
# min_mt_cov:${params.min_mt_cov}
# min_mtcn:${params.min_mtcn}
# max_mtcn:${params.max_mtcn}
# max_contam:${params.max_contam}
# vaf_filter_threshold:${params.vaf_filter_threshold}
# ====================================================
Sample_ID\tStatus\tPass_Filter\tmt_Cov\tmtCN\tContamination\tNote
""".stripIndent()

    """
    # 1. Create the log file with header
    echo "${header}" > pipeline_summary_log.txt
    
    # 2. Combine all individual entry files if they exist
    if [ -n "${entry_files}" ]; then
        cat ${entry_files} | sort >> pipeline_summary_log.txt
    fi

    EXPECTED_COUNT=${expected_count}
    ACTUAL_COUNT=\$(ls *.individual_entry.txt | wc -l)

    # 3. Final check and cleanup
    if [ "${status}" == "FAILED" ] || [ "\$ACTUAL_COUNT" -lt "\$EXPECTED_COUNT" ]; then
        echo -e "\n[ERROR] THE MERGE STAGE FAILED OR SAMPLES MISSING." >> pipeline_summary_log.txt
        echo -e "[ERROR] Expected: \$EXPECTED_COUNT, Found: \$ACTUAL_COUNT" >> pipeline_summary_log.txt
        echo "Error: Stage 2 refinement failed." >&2
        exit 1
    else
        echo -e "\n[INFO] Merging stage completed successfully." >> pipeline_summary_log.txt
        
        # --- NEW CLEANUP LOGIC ---
        # Delete individual entry files only after successful merging
        #echo "[INFO] Cleaning up individual sample files..."
        #rm -f ${entry_files}
    fi
    """
}

// ======================================================================
// WORKFLOW
// ======================================================================
workflow FINALIZER {
  ANNOTATE_COVERAGE(ch_cov_tsvs, ch_hail_dir)
  COMBINE_VCFS(ch_vcfs, ANNOTATE_COVERAGE.out.coverage_mt, ch_hail_dir)
  REFINE_ANNOTATION(COMBINE_VCFS.out.combined_vcf, ch_hail_dir)

  ch_expected_count = ch_vcfs.map { list -> list.size()}
  ch_combine_status = REFINE_ANNOTATION.out.refined_results
        .map { "SUCCESS" }
        .ifEmpty("FAILED")

  COMBINE_SUMMARY_LOGS(ch_summary_entries, ch_combine_status, ch_expected_count)
}

workflow { FINALIZER() }

workflow.onComplete {
  if (workflow.success) {
    log.info "Finalizer completed."
    log.info "Combined MT:   ${params.outdir}/COMBINED_MT/combined_final.mt"
    log.info "Combined VCF:  ${params.outdir}/COMBINED_MT/combined_final.vcf.bgz"
    log.info "Annotations:   ${params.outdir}/combined_annotation"
  } else {
    log.error "Finalizer FAILED. See .nextflow.log and work dirs under: ${workflow.workDir}"
  }
}