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

def vcf_list_raw = all_files.findAll { it.name ==~ /.*\.final\.split\.vcf(\.gz)?$/ }
                            .collect { it.toAbsolutePath().toString() }

def cov_list_raw = all_files.findAll { it.name ==~ /.*\.per_base_coverage\.tsv(\.gz)?$/ }
                            .collect { it.toAbsolutePath().toString() }

if ( vcf_list_raw.isEmpty() ) error "No VCF files found under ${params.merged_dir}"
if ( cov_list_raw.isEmpty() ) error "No Coverage files found under ${params.merged_dir}"

ch_vcf_list_file = Channel.fromList(vcf_list_raw).collectFile(name: 'vcf_paths.txt', newLine: true)
ch_cov_list_file = Channel.fromList(cov_list_raw).collectFile(name: 'cov_paths.txt', newLine: true)

ch_vcf_all = Channel.value(vcf_list_raw)
ch_cov_all = Channel.value(cov_list_raw)

ch_summary_entries = Channel.fromPath("${params.outdir}/*/*.individual_entry.txt").collect()
ch_hail_dir = Channel.fromPath(params.hail_script)
ch_expected_count = Channel.value(vcf_list_raw.size())

// ======================================================================
// PROCESS 1: Convert per-base coverage TSV -> Hail Table/Matrix
// ======================================================================
process ANNOTATE_COVERAGE {
    tag "AnnotateCoverage_MT"
    label 'hail_related'
    publishDir "${params.outdir}/merged_results/coverage_mt", mode: 'copy'

    input:
    path cov_list_file       
    path hail_dir

    output:
    path "coverage.ht", emit: coverage_ht

    script:
    """
    set -euo pipefail

    awk '{
        n = split(\$1, a, "/"); 
        filename = a[n]; 
        split(filename, b, ".");
        sample_id = b[1];
        print sample_id "\t" \$1
    }' ${cov_list_file} > formatted_coverage.list

    ${params.hail_python} ${hail_dir}/annotate_coverage.py \\
      -i formatted_coverage.list \\
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
  path vcf_list_file
  path coverage_ht
  path hail_dir

  output:
  path "combined_final.mt",     emit: combined_mt
  path "combined_final.vcf.bgz", emit: combined_vcf
  path "combined_final.vcf.bgz.tbi", optional: true

  script:
  """
  set -euo pipefail

  awk '{
        n = split(\$1, a, "/"); 
        filename = a[n]; 
        gsub(/\\.merged\\.final\\.split\\.vcf(\\.gz)?\$/, "", filename); 
        print filename "\t" \$1
    }' ${vcf_list_file} > formatted_vcf.list

    ${params.hail_python} ${hail_dir}/combine_vcfs.py \\
    -i formatted_vcf.list \\
    -c ${coverage_ht} \\
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

  input:
  path combined_vcf           
  path pre_list_file         
  path hail_dir              

  output:
  path "annotation", emit: refined_results

  script:
  """
  set -euo pipefail

  mkdir -p tmp_prefiltering
  mkdir -p annotation/final_outputs
  
  CUR_DIR=\$(pwd)
  OUT_DIR="\$CUR_DIR/annotation/final_outputs"
  
  echo "[*] Injecting Sample_ID via path list..."
  tr -d '\\r' < ${pre_list_file} > cleaned_list.txt
  
  while read -r REAL_PATH; do
      [ -z "\$REAL_PATH" ] && continue
      if [ -f "\$REAL_PATH" ]; then
          SAMPLE_NAME=\$(echo \$REAL_PATH | sed "s|${params.outdir}/||" | cut -d'/' -f1)
  
          awk -v sn="\$SAMPLE_NAME" 'BEGIN{FS=OFS="\\t"} {if(NR==1) print "Sample_ID", \$0; else print sn, \$0}' "\$REAL_PATH" > tmp_prefiltering/\${SAMPLE_NAME}.pre.txt
      fi
  done < cleaned_list.txt
  
  python ${hail_dir}/merge_and_calculate.py \\
    --input-dir tmp_prefiltering \\
    --out-dir "\$OUT_DIR" \\
    --pipeline-mode "${params.pipeline_mode}" \\
    --do_post_filter "${params.do_post_filter}"
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
    
    printf "%s\\n" ${entry_files} > file_list.tmp
    
    if [ -s file_list.tmp ]; then
        xargs -a file_list.tmp cat | sort >> pipeline_summary_log.txt
    fi

    EXPECTED_COUNT=${expected_count}
    ACTUAL_COUNT=\$(cat file_list.tmp | wc -l)

    # 4. Final check and cleanup
    if [ "${status}" == "FAILED" ] || [ "\$ACTUAL_COUNT" -lt "\$EXPECTED_COUNT" ]; then
        echo -e "\\n[ERROR] THE MERGE STAGE FAILED OR SAMPLES MISSING." >> pipeline_summary_log.txt
        echo -e "[ERROR] Expected: \$EXPECTED_COUNT, Found: \$ACTUAL_COUNT" >> pipeline_summary_log.txt
        echo "Error: Stage 2 refinement failed." >&2
        exit 1
    else
        echo -e "\\n[INFO] Merging stage completed successfully." >> pipeline_summary_log.txt
    fi
    """
}

// ======================================================================
// WORKFLOW
// ======================================================================
workflow FINALIZER {
    ANNOTATE_COVERAGE(ch_cov_list_file, ch_hail_dir)
    
    COMBINE_VCFS(ch_vcf_list_file, ANNOTATE_COVERAGE.out.coverage_ht, ch_hail_dir)
    
    ch_per_sample_pre_list = Channel.fromPath("${params.outdir}/*/annotation/final_outputs/variant_list_prefiltering.txt")
        .map { it.toAbsolutePath().toString() } 
        .collectFile(name: 'prefiltering_list.txt', newLine: true, sort: true)

    REFINE_ANNOTATION(COMBINE_VCFS.out.combined_vcf, ch_per_sample_pre_list, ch_hail_dir)

    ch_combine_status = REFINE_ANNOTATION.out.refined_results
    .collect()
    .map { "SUCCESS" }
    .ifEmpty("FAILED")

    COMBINE_SUMMARY_LOGS(
    ch_summary_entries, 
    ch_combine_status, 
    ch_expected_count
)
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