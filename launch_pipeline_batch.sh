#!/bin/bash
#SBATCH --job-name=Mito_Pipeline_Batch           # Name of the job
#SBATCH --cpus-per-task=2                        # Recommended: 2 CPUs as one batch handles multiple samples
#SBATCH --mem=8G                                 # Recommended: 8GB depending on BATCH_SIZE
#SBATCH --time=24:00:00                          # Max walltime
#SBATCH --partition=ycga
#SBATCH --output=log/mito_batch_mgr_%j.log       # Main manager log file

# Exit immediately if a command exits with a non-zero status
set -e

# ==============================================================================
#                         User Configuration
# ==============================================================================

# 1) File and Directory Paths
SCRIPT_BASE_DIR="/home/lt692/project_pi_njl27/lt692/human_mito_calling-main"   
MAIN_NF_SCRIPT="$SCRIPT_BASE_DIR/main.nf"                                      
NEXTFLOW_CONFIG="$SCRIPT_BASE_DIR/nextflow_Bouchet.config"                     

MASTER_SAMPLE_LIST="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/test_data/test_sample_cram.tsv"  
OUTPUT_DIR="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/test_data/test_cram_results"             
WORK_DIR_BASE="/home/lt692/scratch_pi_njl27/lt692/human_mt_variant_calling/nextflow_work_cram"                   

# 2) Parallelism and Batch Control
BATCH_SIZE=10            # Number of samples processed within a single Slurm job
CONCURRENT_BATCHES=2     # Maximum number of Batch Jobs allowed to run at once (Slurm throttling)
CLEANUP_ON_SUCCESS=true  # Delete intermediate files after successful completion

# 3) Pipeline Mode Configuration
PIPELINE_MODE="population"
DISEASE_PED_FILE=""

# 4) Environment Setup
module load Nextflow/24.10.2

# Temporary file to store the list of unique samples
INTERNAL_SAMPLE_LIST="${OUTPUT_DIR}/unique_samples_for_job_array.list"

# ==============================================================================
#                            Mode 3: Finalizer Mode (Stage 2)
# ------------------------------------------------------------------------------
# Executed when the script is called with the --finalize argument.
# Merges results from all batches and performs final annotation.
# ==============================================================================
if [ "$1" == "--finalize" ]; then
    echo "========================================================"
    echo "--- STAGE 2: FINALIZER MODE - Generating Combined Results ---"
    echo "========================================================"

    MERGED_INPUT_DIR="${OUTPUT_DIR}/merged_results/wdl_output"
    FINALIZER_WORK_DIR="${WORK_DIR_BASE}/final"
    PROJECT_ROOT=$(dirname "$MAIN_NF_SCRIPT")
    MERGE_NF_SCRIPT="$PROJECT_ROOT/merge.nf"

    # Prepare arguments for disease mode
    EXTRA_ARGS=()
    if [ "${PIPELINE_MODE}" = "disease" ]; then
        EXTRA_ARGS+=( --disease_ped_file "${DISEASE_PED_FILE}" )
    fi

    # Run Nextflow Finalizer
    nextflow run "$MERGE_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        --merged_dir "$MERGED_INPUT_DIR" \
        --input "${INTERNAL_SAMPLE_LIST}" \
        --pipeline_mode "$PIPELINE_MODE" \
        --outdir "$OUTPUT_DIR" \
        "${EXTRA_ARGS[@]}" \
        -w "$FINALIZER_WORK_DIR"

    NF_EXIT=$?

    # Verification: Check if the final variant table was created and is not empty
    FINAL_TABLE=$(find "${OUTPUT_DIR}/merged_results" -name "*variant_list.txt" -size +0c | head -n 1)

    if [ $NF_EXIT -eq 0 ] && [ -n "$FINAL_TABLE" ]; then
        echo "[SUCCESS] Finalizer completed. Output verified: $FINAL_TABLE"
        
        # Cleanup Logic
        if [ "$CLEANUP_ON_SUCCESS" = true ]; then
            echo "[CLEANUP] Removing intermediate and temporary files..."
            [ -d "$WORK_DIR_BASE" ] && rm -rf "$WORK_DIR_BASE"
            find "${OUTPUT_DIR}" -mindepth 1 -type d -name "inputs" -exec rm -rf {} +
            
            # Remove high-volume raw pipeline outputs
            for pipeline_dir in "${OUTPUT_DIR}"/*/variant_calling/MitochondriaMultiSamplePipeline; do
                [ -d "$pipeline_dir" ] && rm -rf "$pipeline_dir"
            done
            echo "[CLEANUP] Finished."
        fi
    else
        echo "[ERROR] Finalizer failed or output missing. Keeping work files for debugging."
        exit 1
    fi
    exit 0

# ==============================================================================
#                            Mode 2: Worker Mode (Stage 1)
# ------------------------------------------------------------------------------
# Executed as part of a Slurm Job Array. 
# Processes a group (Batch) of samples in one Nextflow execution.
# ==============================================================================
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "--- STAGE 1: BATCH WORKER MODE (Batch ID: $SLURM_ARRAY_TASK_ID) ---"

    # Calculate line range for this specific batch
    START_LINE=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE + 1 ))
    END_LINE=$(( START_LINE + BATCH_SIZE - 1 ))

    # Extract sample IDs for this batch from the master list
    CURRENT_BATCH_SAMPLES=$(sed -n "${START_LINE},${END_LINE}p" "$INTERNAL_SAMPLE_LIST")
    
    if [ -z "$CURRENT_BATCH_SAMPLES" ]; then
        echo "[INFO] No samples assigned to this batch index. Exiting."
        exit 0
    fi

    # Define unique work directory for this batch
    BATCH_WORK_DIR="${WORK_DIR_BASE}/batch_job_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$BATCH_WORK_DIR"
    cd "$BATCH_WORK_DIR"

    # Create a batch-specific TSV input file for Nextflow
    BATCH_TSV="${BATCH_WORK_DIR}/batch_${SLURM_ARRAY_TASK_ID}_input.tsv"
    > "$BATCH_TSV" # Initialize empty file
    
    for SID in $CURRENT_BATCH_SAMPLES; do
        awk -v id="$SID" '$1==id {print}' "$MASTER_SAMPLE_LIST" >> "$BATCH_TSV"
    done

    echo "[INFO] Processing the following samples in this batch:"
    echo "$CURRENT_BATCH_SAMPLES"

    # Run Nextflow on the batch TSV
    nextflow run "$MAIN_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        -w "$BATCH_WORK_DIR" \
        --input "$BATCH_TSV" \
        --outdir "$OUTPUT_DIR" \
        --pipeline_mode "$PIPELINE_MODE" \
        $( [ "$PIPELINE_MODE" == "disease" ] && echo "--disease_ped_file $DISEASE_PED_FILE" )

    exit 0

# ==============================================================================
#                            Mode 1: Master Mode
# ------------------------------------------------------------------------------
# Main entry point. Calculates batches and submits Stage 1 and Stage 2 jobs.
# ==============================================================================
else
    echo "--- STAGE 0: MASTER MODE - Initializing Batching and Submission ---"
    mkdir -p log
    mkdir -p "$OUTPUT_DIR"
    
    # 1. Generate a unique sample list (excluding header)
    cut -f1 "$MASTER_SAMPLE_LIST" | grep -vi "^sample" | sort -u > "$INTERNAL_SAMPLE_LIST"
    NUM_SAMPLES=$(wc -l < "$INTERNAL_SAMPLE_LIST")

    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "[ERROR] No samples found in $MASTER_SAMPLE_LIST"
        exit 1
    fi

    # 2. Calculate the number of batches required
    # Formula: ceil(total / batch_size)
    NUM_BATCHES=$(( (NUM_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE ))
    MAX_ARRAY_INDEX=$(( NUM_BATCHES - 1 ))

    echo "[INFO] Total Samples Found: $NUM_SAMPLES"
    echo "[INFO] Samples per Job (Batch Size): $BATCH_SIZE"
    echo "[INFO] Total Slurm Jobs to submit: $NUM_BATCHES"
    echo "[INFO] Max Concurrent Running Batches: $CONCURRENT_BATCHES"

    # 3. Submit Stage 1 Job Array 
    # The %$CONCURRENT_BATCHES syntax limits the number of simultaneously running tasks
    STAGE1_JOB_ID=$(sbatch --parsable \
        --array=0-${MAX_ARRAY_INDEX}%${CONCURRENT_BATCHES} \
        --output="log/mito_worker_batch_%A_%a.log" \
        "$0")
    echo "[INFO] Stage 1 (Batch Array) Job ID: ${STAGE1_JOB_ID}"

    # 4. Submit Stage 2 Finalizer 
    # Dependency 'afterok' ensures this only runs after all Stage 1 tasks succeed
    STAGE2_JOB_ID=$(sbatch --parsable \
        --dependency=afterok:${STAGE1_JOB_ID} \
        --output="log/mito_finalizer_%A.log" \
        "$0" --finalize)
        
    echo "[INFO] Stage 2 (Finalizer) Job ID: ${STAGE2_JOB_ID}"
    echo "[DONE] Workflow submitted. Monitor status using 'squeue -u $USER'"
fi