#!/bin/bash
#SBATCH --job-name=Mito_Pipeline_Manager         # Job name in the Slurm queue
#SBATCH --cpus-per-task=2                        # Request 2 CPU cores
#SBATCH --mem=8G                                 # Request 8 GB of memory
#SBATCH --time=24:00:00                          # Max runtime of 24 hours
#SBATCH --output=log/mito_pipeline_manager_%j.log # Log file name (%j = job ID)

# Exit immediately if any command fails
set -e

# ==============================================================================
#                         User Configuration
# ==============================================================================

# 1) File and Directory Paths
SCRIPT_BASE_DIR="/path/to/human_mito_calling-main"   # Pipeline root
MAIN_NF_SCRIPT="$SCRIPT_BASE_DIR/main.nf"                                      # Main Nextflow script
NEXTFLOW_CONFIG="$SCRIPT_BASE_DIR/nextflow_Bouchet.config"                     # Nextflow config

MASTER_SAMPLE_LIST="/path/to/test_sample_cram.tsv"  # TSV with all samples/URLs
OUTPUT_DIR="/path/to/output_dir"             # Final output dir (fallback)
WORK_DIR_BASE="/path/to/work_dir"                   # Base work dir for intermediates                  

# 2) Cleanup Configuration
CLEANUP_ON_SUCCESS=true               # Set to 'true' to delete intermediate work and temporary input directories on success

# 3) Slurm Job Array Concurrency
CONCURRENT_SAMPLES=3                     # Max number of array tasks to run simultaneously

# 4) Pipeline Mode: "population" or "disease"
PIPELINE_MODE="disease"

# Required only if PIPELINE_MODE="disease"; ignored for "population"
DISEASE_META_FILE="/path/to/disease_meta.tsv"

# 5) Environment for Nextflow (please change this per your server)
module load Nextflow/24.04.4

# ==============================================================================
#                                  Script Modes
# ------------------------------------------------------------------------------
# 1) Master Mode (default, no args, login node):
#       - Build unique sample list
#       - Submit Stage 1 (array workers) and Stage 2 (finalizer) jobs
# 2) Worker Mode (compute node; when $SLURM_ARRAY_TASK_ID is set):
#       - Process one unique sample (one array task)
# 3) Finalizer Mode (compute node; when $1 == --finalize):
#       - Merge/annotate final outputs across samples (Hail/Nextflow)
# ==============================================================================

# Internal file (do not change)
INTERNAL_SAMPLE_LIST="unique_samples_for_job_array.list"

# ==============================================================================
#                            Mode 3: Finalizer Mode (Stage 2)
# ==============================================================================
if [ "$1" == "--finalize" ]; then
    echo "========================================================"
    echo "--- STAGE 2: FINALIZER MODE - Generating Combined Results ---"
    echo "========================================================"

    MERGED_INPUT_DIR="${OUTPUT_DIR}/merged_results/wdl_output"
    FINALIZER_WORK_DIR="${WORK_DIR_BASE}/final"
    PROJECT_ROOT=$(dirname "$MAIN_NF_SCRIPT")
    MERGE_NF_SCRIPT="$PROJECT_ROOT/merge.nf"

    # Assemble disease mode arguments
    EXTRA_ARGS=()
    if [ "${PIPELINE_MODE}" = "disease" ]; then
        EXTRA_ARGS+=( --disease_meta_file "${DISEASE_META_FILE}" )
    fi

    # Run Nextflow Finalizer
    nextflow run "$MERGE_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        --merged_dir "$MERGED_INPUT_DIR" \
        --pipeline_mode "$PIPELINE_MODE" \
        --outdir "$OUTPUT_DIR" \
        "${EXTRA_ARGS[@]}" \
        -w "$FINALIZER_WORK_DIR"

    NF_EXIT=$?

    # --- Physical Validation: Check if the final variant table exists and is not empty ---
    FINAL_TABLE=$(find "${OUTPUT_DIR}/merged_results" -name "*variant_list.txt" -size +0c | head -n 1)

    if [ $NF_EXIT -eq 0 ] && [ -n "$FINAL_TABLE" ]; then
        echo "[SUCCESS] Finalizer completed successfully. Output verified: $FINAL_TABLE"
        
        # --- Cleanup Logic ---
        if [ "$CLEANUP_ON_SUCCESS" = true ]; then
            echo "[CLEANUP] Cleaning up intermediate and temporary files..."
            
            # 1. Remove the base Nextflow work directory (Intermediate computation files)
            if [ -d "$WORK_DIR_BASE" ]; then
                echo "[CLEANUP] Removing work directory: $WORK_DIR_BASE"
                rm -rf "$WORK_DIR_BASE"
            fi
            
            # 2. Remove all 'inputs' directories within the output directory (Temporary TSVs)
            find "${OUTPUT_DIR}" -mindepth 1 -type d -name "inputs" -exec echo "[CLEANUP] Removing: {}" \; -exec rm -rf {} +

            # 3. Remove individual sample entry text files from results
            echo "[CLEANUP] Removing individual entry logs from output..."
            find "${OUTPUT_DIR}" -name "*.individual_entry.txt" -type f -delete

            # 4. Remove the raw WDL output folder to save significant space
            # Path: [OUTPUT_DIR]/variant_calling/MitochondriaMultiSamplePipeline
            echo "[CLEANUP] Searching for raw pipeline outputs in: ${OUTPUT_DIR}/*/variant_calling/MitochondriaMultiSamplePipeline"
            
            # Use a for loop to handle the wildcard expansion safely
            for pipeline_dir in "${OUTPUT_DIR}"/*/variant_calling/MitochondriaMultiSamplePipeline; do
                if [ -d "$pipeline_dir" ]; then
                    echo "[CLEANUP] Removing raw pipeline output: $pipeline_dir"
                    rm -rf "$pipeline_dir"
                fi
            done
            echo "[CLEANUP] Cleanup finished."
        fi
    else
        echo "[ERROR] Finalizer failed or critical output missing. Keeping work files for debugging."
        exit 1
    fi
    exit 0

# ==============================================================================
#                            Mode 2: Worker Mode (Stage 1)
# ==============================================================================
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "--- STAGE 1: WORKER MODE - Processing Unique Sample ---"

    TASK_INDEX=$((SLURM_ARRAY_TASK_ID + 1))
    SAMPLE_ID=$(awk -v line=$TASK_INDEX 'NR==line {print; exit}' "$INTERNAL_SAMPLE_LIST" | tr -d '\r')

    WORK_DIR="${WORK_DIR_BASE}/batch_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    TEMP_TSV="${WORK_DIR}/sample_${SAMPLE_ID}.tsv"
    awk -v id="$SAMPLE_ID" '$1==id {print}' "$MASTER_SAMPLE_LIST" > "$TEMP_TSV"

    nextflow run "$MAIN_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        -w "$WORK_DIR" \
        --input "$TEMP_TSV" \
        --outdir "$OUTPUT_DIR" \
        --pipeline_mode "$PIPELINE_MODE" \
        $( [ "$PIPELINE_MODE" == "disease" ] && echo "--disease_meta_file $DISEASE_META_FILE" )

    exit 0

# ==============================================================================
#                            Mode 1: Master Mode
# ==============================================================================
else
    echo "--- STAGE 0: MASTER MODE - Submitting Job Array ---"
    mkdir -p log
    cut -f1 "$MASTER_SAMPLE_LIST" | grep -vi "^sample" | sort -u > "$INTERNAL_SAMPLE_LIST"
    NUM_SAMPLES=$(wc -l < "$INTERNAL_SAMPLE_LIST")
    ARRAY_INDEX=$((NUM_SAMPLES - 1))

    STAGE1_JOB_ID=$(sbatch --parsable --array=0-${ARRAY_INDEX}%${CONCURRENT_SAMPLES} \
        --output="log/mito_pipeline_worker_%A_%a.log" "$0")
    echo "[INFO] Stage 1 Job Array ID: ${STAGE1_JOB_ID}"

    STAGE2_JOB_ID=$(sbatch --parsable --dependency=afterok:${STAGE1_JOB_ID} \
        --output="log/mito_pipeline_finalizer_%A.log" "$0" --finalize)
    echo "[INFO] Stage 2 Finalizer Job ID: ${STAGE2_JOB_ID}"
fi