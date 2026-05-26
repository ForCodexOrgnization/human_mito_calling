# Human mtDNA Variant Calling Pipeline
This Nextflow pipeline is designed to perform human mitochondrial (mtDNA) variant calling, starting from raw FASTQ/BAM/CRAM files. It utilizes the GATK mutect2 mitochondrial-mode workflow (implemented in WDL and executed via Cromwell) and a subsequent Hail-based annotation workflow.

The pipeline supports both "population" and "disease" analysis modes, allowing for flexible handling of different study designs and objectives.

The pipeline is primarily configured to run on a Slurm cluster environment using a three-stage approach managed by an SLURM Job Array Manager.

## Quick Start
1. Clone the repository:
```bash
git clone https://github.com/LiTang1013/human_mito_calling.git
cd human_mito_calling
``` 
2. Make sure you have all prerequisites installed (see [Prerequisites](#prerequisites)).
3. Prepare your sample list TSV file (e.g., `test_sample.tsv`), following the specified format (see [Input Files](#input-files)). If you want to use the "disease" mode, prepare the disease metadata file (e.g. `disease_meta.tsv`) as well.
4. Prepare required files for the pipeline (see [Required Files](#required-files)).
5. Review and update paths in `nextflow.config` to match your environment. As well as change sample-level/variant filtering parameters if needed (see [Configuration](#configuration)).
6. Review and update paths in `launch_pipeline.sh` and set launch parameters (see [Launch](#launch)).
7. Launch the pipeline using the provided Slurm submission script, make sure you are in the base directory of the repository:
```bash
bash launch_pipeline.sh
```
8. Monitor the Slurm job logs for progress and any potential issues (see [Monitoring](#monitoring)).
9. After a successful run, review the output directory structure (see [Output Structure](#output-structure)).

### Note
Processing three CRAM test samples (~5.3 GB each) completes in about 48m 57s end-to-end with the default settings (submit-rate limit 180/h, queue size 20; see [Configuration](#configuration) for CPU/memory). You can adjust these parameters to match your cluster capacity.


## Prerequisites
Before running the pipeline, ensure you have the following installed and configured:
* Nextflow (v24.10.2 or later recommended)
* Cromwell (You can download it [here](https://github.com/broadinstitute/cromwell/releases/download/90/cromwell-90.jar); only the JAR file is needed)
* haplocheck (You can download it [here](https://github.com/leklab/haplocheckCLI?tab=readme-ov-file); only the JAR file is needed)
* Conda/Miniconda (You will need a dedicated conda environment. We recommend using our [install script](https://github.com/LiTang1013/Human_mito_calling/edit/main/install_env.sbatch) to create the `bio_pipeline_env`)
    Within this environment, the following software should be installed:
    - bwa
    - samtools
    - picard
    - gatk4
    - python (with hail, pandas and numpy packages)
    - R
    - openjdk
    - VEP

* Slurm Workload Manager (required for the cluster profile)


## Input Files
The pipeline requires a sample list TSV file (e.g., test_sample.tsv) and a disease metadata file (e.g. `disease_meta.tsv`), which is only needed if you want to use the "disease" mode. The sample list file could contain multiple samples with FASTQ/BAM/CRAM paths or URLs. Each sample can have multiple sequencing runs, the pipeline will merge them as one in the processing steps (if you don't want to merge them, you can specify them with different sample IDs).

### Sample list TSV file
A tab-separated file with three columns.

- For FASTQ input:
```
SampleID	R1_Path	R2_Path
Sample1	/path/to/data/Sample1_L1_R1.fastq.gz	/path/to/data/Sample1_L1_R2.fastq.gz
Sample1	/path/to/data/Sample1_L2_R1.fastq.gz	/path/to/data/Sample1_L2_R2.fastq.gz
Sample2	/path/to/data/Sample2_L1_R1.fastq.gz	/path/to/data/Sample2_L1_R2.fastq.gz
Sample2	/path/to/data/Sample2_L2_R1.fastq.gz	/path/to/data/Sample2_L2_R2.fastq.gz
```
- For BAM input:
```
SampleID	BAM_Path	Index_Path
Sample1	/path/to/data/Sample1.bam   /path/to/data/Sample1.bai
Sample2	/path/to/data/Sample2.bam   /path/to/data/Sample2.bai
```
- For CRAM input:
```
SampleID	CRAM_Path	Index_Path
Sample1	/path/to/data/Sample1.cram   /path/to/data/Sample1.crai
Sample2	/path/to/data/Sample2.cram   /path/to/data/Sample2.crai
``` 

### Disease metadata file format
Only needed for "disease" mode. A tab-separated file with at least the following columns.
```
SampleID	FamilyID	Category	OtherInfo
Sample1	Family1	proband	SomeInfo1
Sample2	Family1	control	SomeInfo2
```
 Must contain at least the following columns: `SampleID`, `FamilyID`, `Category` (proband/control), and any other metadata columns you wish to include.

## Required Files
- Whole genome reference and indices

If you are using FASTQ input, the pipeline requires a complete reference genome FASTA file with its associated index and dictionary files. If you are using BAM/CRAM input, ensure that the files are aligned to the same reference genome specified in `nextflow.config`. The whole genome reference files will be used for alignment and mtCN estimation.

- Mitochondrial reference and indices

This includes the standard mitochondrial reference FASTA and a shifted version (shifted by 8000 bases) along with their BWA indices. These are necessary for the GATK mutect2 mitochondrial-mode workflow.

- Other required files for variant calling and filtering

Additional files such as the shift back chain file and blacklisted sites BED files are also required for accurate variant calling and filtering.

`autosome_XY.interval_list` is used to calculate the autosomal coverage for mtCN estimation.

Here we provide an example set of reference and other required files needed for the pipeline (hg38/GRCh38) in the `required_files` directory:
```
/Homo_sapiens_assembly38.fasta
/Homo_sapiens_assembly38.fasta.fai
/Homo_sapiens_assembly38.dict
/Homo_sapiens_assembly38.chrM
/Homo_sapiens_assembly38.chrM.fai
/Homo_sapiens_assembly38.chrM.dict
/Homo_sapiens_assembly38.chrM.fasta.amb
/Homo_sapiens_assembly38.chrM.fasta.bwt
/Homo_sapiens_assembly38.chrM.fasta.pac
/Homo_sapiens_assembly38.chrM.fasta.sa
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.dict
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac
/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa
ShiftBack.chain
blacklist_sites.hg38.chrM.bed
non_control_region.interval_list
control_region_shifted.interval_list
autosome_XY.interval_list
```
- Annotation database files

Here we use some common annotation databases for the pipeline to annotate mtDNA variants, such as ClinVar, gnomAD, HelixMTdb, MITOMAP, MitoTIP and so on. We provide the complete annotation database files for [download](https://example.com/download). The paths of these files should be specified in the `nextflow.config`. Note that the human VEP cache files are based on v113, if your VEP version is different, you may need to update the cache files accordingly (the appropriate version can be downloaded from the [VEP FTP](https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache)).


## Configuration

The pipeline is configured via `nextflow.config`. This file controls Nextflow's software paths, and resource defaults.
<div style="overflow-x: auto;">

| **Parameter Group** | **Key** | **Description**  | **User Edit?** |
|----------------------|----------|-----------------|----------------|
| **Conda Settings** | `conda_base_path` | Path to your Conda base installation <br> *(Run `conda info` to check the `base environment` path)* | Yes |
|  | `conda_env_path` | Path to the Conda environment used by the pipeline. <br> *(Run `conda info` to check the `active env location` path)*  | Yes |
| **Core Software**  | `cromwell_jar` | Path to the Cromwell JAR file. | Yes |
|  | `software_base_dir` | Base directory for haplocheck and cromwell. | Yes |
| **Required Files** | `required_base_dir` | Base directory for all required reference genomes, indexes, and annotation files.| Yes |
| **Filtering Parameters** | `vaf_filter_threshold` | The minimum filter threshold of allele frequency for variant calling. | Optional <br> *(Default: `0.01`)* |
|   | `min_mt_cov` | Minimum mitochondrial coverage for sample-level filtering. | Optional <br> *(Default: `100`)* |
|   | `min_mtcn` | Minimum mitochondrial copy number for sample-level filtering. | Optional <br> *(Default: `50`)* |
|   | `max_mtcn` | Maximum mitochondrial copy number for sample-level filtering. | Optional <br> *(Default: `500`)* |
|   | `max_contam` | Maximum contamination level (%) for sample-level filtering. | Optional <br> *(Default: `0.02`)* |
| **WDL Config** | `params.wdl_inputs` | Pass the required software and file paths to the WDL. *(If your files are already in the `required_files` folder, no changes are needed)*| Optional |
| **Hail Config** | `params.hail_pipeline_config` | A map for Hail-related settings, including VEP cache directory, annotation databases, and metadata files. *(If your files are already in the `annotation_databases` folder, no changes are needed)*| Optional |
| **Execution Profiles** |`params.cromwell_options`| Parameters controlling how Cromwell submits its own jobs to Slurm.| Optional *(Change according to your Slurm setup)* |
| | `process.withLable`| Resource allocations for each process.| Optional *(Change according to your Slurm setup)* |
| | `profiles`| Uses the Slurm executor and enables Conda. <br> Automatically activates the environment before each job using and sets the internal WDL backend to Slurm| No |
| **Executor Config** | `executor` | Limits the job submission rate to avoid overwhelming the scheduler.| Optional *(Change according to your Slurm setup)* |

</div>


## Launch
The entire pipeline is launched using the provided Slurm submission script. The `launch_pipeline.sh` script will perform the following automatically:

- Master Mode: Calculate the total number of samples and submit the Stage 1 Job Array.

- Worker Mode (Stage 1): Each array task will run the Nextflow pipeline in a dedicated work directory for one sample to process the preprocessing, variant calling, and initial annotation steps.

- Finalizer Mode (Stage 2): A final job is submitted with an afterok dependency on the array, executing the Hail merging and final annotation only after all sample-specific jobs succeed.

Before launching, review and update the following variables in `launch_pipeline.sh` to match your environment:

<div style="overflow-x: auto;">

| Variable | Description | Default Path (Needs Review) |
|-----------|--------------|------------------------------|
| SCRIPT_BASE_DIR | Base directory of the pipeline repository/scripts | `/path/to/human_mito_calling` |
| MAIN_NF_SCRIPT | Stage 1 Nextflow pipeline script (*main.nf*). | `/path/to/main.nf` |
| FINALIZE_NF_SCRIPT | Stage 2 finalizer workflow (*merge.nf*). | `/path/to/merge.nf` |
| NEXTFLOW_CONFIG | Nextflow configuration file (*nextflow.config*). | `/path/to/nextflow.config` |
| MASTER_SAMPLE_LIST | Master sample list TSV | `/path/to/test_sample_list.tsv` |
| OUTPUT_DIR | Root output directory (Stage 1 & Stage 2 results). | `/path/to/results` |
| WORK_DIR_BASE | Base working directory for all jobs. | `/path/to/workdir` |
| CONCURRENT_SAMPLES | Number of concurrent samples to process in the Slurm Job Array. | `3` |
| PIPELINE_MODE | Pipeline analysis mode: "population" or "disease". | `population` |
| DISEASE_META_FILE | Disease metadata TSV file (only needed for "disease" mode). | `/path/to/disease_meta.tsv` |
| Module | Load required module (Or make sure your Nextflow environments is set up) | `module load Nextflow/24.10.2` |

</div>

## Pipeline Overview
Will be added later (maybe a figure).


## Monitoring
You can monitor the progress of the pipeline by checking the Slurm job status using the `squeue` command. Additionally, log files for each job will be generated in the `log` folder with the same directory as the script (e.g., `/log/mito_pipeline_manager_xxxxx_0.log` for each job). The log files record job submissions, execution status, and any errors encountered. An example of a successful run log is shown below:
```
=================================================================
--- STAGE 1: WORKER MODE on Compute Node (Task 0) ---
=================================================================
[*] Processing unique sample: SampleID
[*] Changed to unique work/run directory: /path/to/workdir/batch_0/
[*] Created/Updated stable input file: /path/to/workdir/batch_0/sample_xxx.tsv
[*] Input file contents:
SampleID	R1_Path	R2_Path
[*] Running in 'disease' mode. Checking for Disease meta file.
[*] Starting Nextflow...

PIPELINE START
=================================
Input file:        /path/to/workdir/batch_0/sample_xxx.tsv
Output directory:  /path/to/results
=================================

executor >  slurm (1)
[7f/005e0b] ALI…A-MEM on HG00258 (Pair 8)) | 10 of 10, cached: 10 ✔
[00/1f21f9] SOR…CRAM for HG00258 (Pair 8)) | 10 of 10, cached: 10 ✔
[53/7d6cbc] CONVERT_BAM_TO_CRAM            -
[a4/6d3041] MER… (Merge CRAMs for HG00258) | 1 of 1, cached: 1 ✔
[94/6fc315] GEN…TSV (CRAM TSV for HG00258) | 1 of 1, cached: 1 ✔
[c6/d103f0] GEN…SON (WDL JSON for HG00258) | 1 of 1, cached: 1 ✔
[3a/6abd17] RUN… (Variant Calling HG00258) | 1 of 1, cached: 1 ✔
[89/7d6cbc] CAL…alculate_mtCN for HG00258) | 1 of 1, cached: 1 ✔
[fd/7d6cbc] SAMPLE_LEVEL_FILTER            | 1 of 1, cached: 1 ✔
[s3/7d6cbc] ANNOTATE_INDIVIDUAL_VCF        | 1 of 1, cached: 1 ✔

Pipeline completed successfully.
Output files are in: /path/to/result

Completed at: 10-Nov-2025 13:16:00
Duration    : 1h 22m
CPU hours   : 18.9 (98.5% cached)
Succeeded   : 1
Cached      : 12

Batch 0 completed successfully for sample xxx.
--- Finished Job Array Task 0 ---
```
Each sample-specific job will have its own work directory (e.g., `/path/to/workdir/batch_0/`), where you can find detailed Nextflow logs and intermediate files for troubleshooting if needed. 
 The structure of the work directory will look like this:
```
/path/to/workdir/
├── batch_0/            # Sample-specific work directory
│   ├── 7f/            # Directories for each Nextflow process
│   │   ├── 005e0b85dd90454ccc108d80e9800a/       # Temporary directory for intermediate files
│   │   └── .command.log                          # Detailed log file for the process
│   └── sample_xxx.tsv                            # Stable input file for the sample
├── batch_1/
│   └── ...
└── ... (other batches)
```
The two-character hash (e.g., 7f) indicates specific process, you can use this hash value to locate the directory of each process. If you want to check the detailed log for a specific process, use following:
```
cat /path/to/workdir/batch_0/7f/005e0b85dd90454ccc108d80e9800a/.command.log
```
Once you’ve fixed the issues, just run `bash ./launch_pipeline.sh`. Thanks to Nextflow’s resume mechanism, the pipeline will pick up exactly where it failed. All intermediate outputs and task statuses are cached in the work/ directory (e.g., `/path/to/workdir/batch_0/`), so successful steps are not re-executed—only the failed or newly affected tasks will run.

### A few tips

- When it resumes: If inputs, params, or scripts for a task haven’t changed, Nextflow reuses the cached result. If any of those change, that task (and its dependents) are re-run automatically.

- Forcing a clean re-run of specific steps: delete the corresponding subfolders under work/ (or bump a version param) and rerun.

- Full fresh run: remove work/ and the output directory (or run with a new -w and --outdir).

## Output Structure
All results are deposited in the directory specified by `OUTPUT_DIR` in the `launch_pipeline.sh`.
The structure is:
```
/path/to/results/
├── Sample1/
│   ├── alignment/                 # Only created when the input is FASTQ
│   │   ├── Sample1.merged.cram
│   │   └── Sample1.merged.cram.crai
│   ├── variant_calling/
│   │   ├── Sample1.final.vcf.gz
│   │   ├── Sample1.per_base_coverage.tsv
│   │   └── ... (other WDL outputs)
│   ├── mtCN/
│   │   └── mtCN_summary.txt
│   └── annotation/
│       ├── vep_vcf/
│       │   └── Sample1.merged_vep.vcf
│       ├── metadata/
│       │   ├── contamination.txt
│       │   └── haplogroup_full.txt
│       └── final_outputs/
│           ├── variant_list.txt
│           └── variant_list_prefiltering.txt
├── ... (Other samples)
│
└── merged_results                  # Stage 2 Hail Final Outputs
    ├── wdl_output/                 # The wdl outputs of all samples required by Hail
    │   ├── Sample1.merged.final.split.vcf
    │   ├── Sample1.merged.haplocheck_contamination.txt
    │   ├── Sample1.merged.per_base_coverage.tsv
    │   └── ...
    ├── coverage_mt/                # Hail Matrix Table for coverage
    │   ├── coverage.ht
    │   ├── coverage.mt
    │   ├── coverage_files.list
    ├── combined_mt/                
    │   ├── combined_final.mt       # Combined Hail Matrix Table
    │   └── combined_final.vcf.bgz  # Final annotated combined VCF
    └── annotation/
        ├── vep_vcf/
        │   └── combined_final.filtered.vcf  # Will be "combined_final.proband_only.filtered.vcf", under "disease" mode
        │   └── combined_final.filtered_vep.vcf # Will be "combined_final.proband_only.filtered_vep.vcf", under "disease" mode
        └── final_outputs/
            ├── variant_list.txt    # Will be "Proband_variant_list.txt", under "disease" mode
            └── variant_list_prefiltering.txt  # Will be "Proband_variant_list_prefiltering.txt", under "disease" mode             

```
