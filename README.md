# RNAseq-pipeline

A Nextflow DSL2 pipeline for end-to-end RNA-seq analysis — from raw FASTQ reads to differential expression results and publication-quality figures.

**Author:** Kalpani De Silva (kdesilva@wi.mit.edu)
**Group:** Bioinformatics and Research Computing (BaRC), Whitehead Institute
**Version:** 2.0-portable — May 2026
**Location:** `/nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline`


---

## Pipeline Overview

```
Raw FASTQ → FastQC → Trim Galore → FastQC (post-trim)
         → STAR alignment
         → infer_experiment.py (auto strandedness on all BAMs, optional)
         → featureCounts
         → merge_counts.py (combine per-sample counts)
         → id_to_name.py (gene ID → gene name, auto-extracted from GTF)
         → RunDESeq2.R (DESeq2 differential expression)
         → cluster_draw_pheatmap.R (heatmap of significant genes)
         → MultiQC (aggregate QC report)
```

| Step | Tool | Description |
|------|------|-------------|
| 1 | FastQC | Raw read quality control |
| 2 | Trim Galore | Adapter trimming and quality filtering |
| 3 | FastQC | Quality control on trimmed reads |
| 4 | STAR | Splice-aware alignment to reference genome |
| 5 | RSeQC `infer_experiment.py` | Auto-detect library strandedness on every BAM *(optional, requires `--bed`)* |
| 6 | featureCounts | Gene-level read quantification |
| 7 | `merge_counts.py` | Merge per-sample counts into a single matrix |
| 8 | `id_to_name.py` | Convert Ensembl gene IDs to gene names *(automatic — extracted from GTF; use `--annotation` to override with a custom mapping file)* |
| 9 | `RunDESeq2.R` | Differential expression with DESeq2 |
| 10 | `cluster_draw_pheatmap.R` | Clustered heatmap of significant DE genes |
| 11 | MultiQC | Aggregate QC report across all samples |

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.10
- Java 11 or later
- Singularity (cluster) **or** Conda (laptop/local)

---

## Software Environment

All tools are bundled in a single self-contained Singularity image (`singularity_cache/rnaseqpipe-portable-2.0.sif`):

| Category | Tools |
|----------|-------|
| QC | FastQC, Trim Galore, MultiQC |
| Alignment | STAR 2.7.1a, samtools |
| Quantification | featureCounts (Subread), RSeQC |
| R / Bioconductor | R 4.2+, DESeq2, vsn, ComplexHeatmap, pheatmap, ggplot2, ggrepel |
| Python | Python 3.10+, pandas |
| Clustering | Cluster 3.0 |

### Singularity container setup

**Whitehead cluster users:** The pre-built image is already present — no build step needed.

**External users — copy the image** (simplest):

```bash
scp <whitehead-user>@<cluster>:/nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/singularity_cache/rnaseqpipe-portable-2.0.sif \
    /path/to/your/RNAseq-pipeline/singularity_cache/
```

**External users — rebuild the image** (if copying isn't possible):

```bash
cd /path/to/your/RNAseq-pipeline
singularity build --fakeroot singularity_cache/rnaseqpipe-portable-2.0.sif \
                             singularity_cache/Singularity.def
```

Use `--fakeroot` on shared clusters where you do not have root access. The build takes 20–40 minutes and produces a ~2–4 GB `.sif` file.

---

## Quick Start

### External users (any system with Singularity)

Copy or rebuild the `.sif` (see above), then run with `-profile singularity`:

```bash
nextflow run /path/to/RNAseq-pipeline/main.nf \
  -profile singularity \
  --reads        "/path/to/fastq/*_{R1,R2}.fastq.gz" \
  --gtf          /path/to/annotation.gtf \
  --star_index   /path/to/star_index \
  --samplesheet  /path/to/samplesheet.csv \
  --control_group WT \
  --bed          /path/to/annotation.bed \
  --outdir results
```

On an external HPC with SLURM, use `-profile cluster` instead and update `conf/slurm.config` with your cluster's partition name and account (see Execution Profiles below).

### Minimal run — QC + alignment only, no DESeq2 (Whitehead cluster)

```bash
nextflow run /nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/main.nf \
  --reads      "/path/to/data/*_{R1,R2}.fastq.gz" \
  --gtf        /path/to/annotation.gtf \
  --star_index /path/to/star_index \
  --skip_deseq2 \
  --outdir results \
  -profile cluster
```

### Full run with differential expression (Whitehead cluster)

```bash
nextflow run /nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/main.nf \
  --reads        "/path/to/data/*_{R1,R2}.fastq.gz" \
  --gtf          /path/to/annotation.gtf \
  --star_index   /path/to/star_index \
  --samplesheet  /path/to/samplesheet.csv \
  --control_group WT \
  --bed          /path/to/annotation.bed \
  --outdir results \
  -profile cluster
```

### Override gene name mapping (optional)

Gene names are extracted automatically from the GTF. To use a custom mapping file instead (e.g. from Ensembl BioMart):

```bash
nextflow run /nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/main.nf \
  --reads          "/path/to/data/*_{R1,R2}.fastq.gz" \
  --gtf            /path/to/annotation.gtf \
  --star_index     /path/to/star_index \
  --samplesheet    /path/to/samplesheet.csv \
  --control_group  WT \
  --bed            /path/to/annotation.bed \
  --annotation     /path/to/gene_id_to_name.txt \
  --annotation_key "Gene name" \
  --outdir results \
  -profile cluster
```

### Single-end reads

```bash
nextflow run /nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/main.nf \
  --reads      "/path/to/data/*.fastq.gz" \
  --single_end \
  --gtf        /path/to/annotation.gtf \
  --star_index /path/to/star_index \
  --samplesheet samplesheet.csv \
  --outdir results \
  -profile cluster
```

### Resume a failed or interrupted run

Add `-resume` to any command to reuse cached results from completed steps:

```bash
nextflow run ... -resume
```

---

## Execution Profiles

| Profile | Executor | Container | Use case |
|---------|----------|-----------|----------|
| `standard` | local | none | Quick local test (tools must be in PATH) |
| `conda` | local | conda env | Local run via conda (R + Bioconductor included) |
| `singularity` | local | `.sif` | Laptop or external cluster with Singularity |
| `cluster` | SLURM | `.sif` | Whitehead HPC (partition `24`, Singularity); update `conf/slurm.config` for other clusters |

### Running on the Whitehead HPC (SLURM)

The `cluster` profile submits all jobs to SLURM partition `24` under account `wibrusers`, with Singularity as the container runtime. Run the Nextflow head process itself on a compute node to avoid running on the head node:

```bash
srun --partition=24 --account=wibrusers --cpus-per-task=2 --mem=8G --time=24:00:00 \
  nextflow run /nfs/BaRC_Public/BaRC_code/pipelines/RNAseq-pipeline/main.nf \
    --reads      "/path/to/fastq/*_{R1,R2}.fastq.gz" \
    --gtf        /path/to/annotation.gtf \
    --star_index /path/to/star_index \
    --samplesheet samplesheet.csv \
    --outdir results \
    -profile cluster \
    -resume
```

---

## Parameters

### Input / Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `*.fastq.gz` | Glob pattern for input FASTQ files (must be quoted) |
| `--single_end` | `false` | Set to `true` for single-end reads |
| `--outdir` | `results` | Output directory |
| `--gtf` | *(required)* | Genome annotation GTF |
| `--star_index` | *(required)* | Pre-built STAR index directory |

### Strandedness

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bed` | `null` | BED file for `infer_experiment.py`. If provided, strandedness is auto-detected across all BAMs (majority vote). |
| `--fc_strandedness` | `0` | Manual strandedness (only when `--bed` is not set): `0`=unstranded, `1`=forward, `2`=reverse |

### Trimming

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_length` | `20` | Minimum read length after trimming |
| `--quality_cutoff` | `20` | Phred quality cutoff |
| `--skip_trimming` | `false` | Skip Trim Galore |

### Alignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `8` | Threads for STAR / featureCounts |
| `--star_ram` | `40` | Max RAM (GB) for STAR |
| `--star_twopass` | `false` | Enable STAR 2-pass mode |

### Differential Expression

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `null` | CSV with sample metadata (see format below) |
| `--condition_col` | `condition` | Column in samplesheet for DE contrast |
| `--control_group` | `null` | Reference/control condition name (e.g. `WT`). If not set, the first condition in the samplesheet is used. |
| `--lfc_threshold` | `1.0` | \|log2FC\| threshold for plots and heatmap filtering |
| `--padj_cutoff` | `0.05` | FDR threshold for plots and heatmap filtering |
| `--skip_deseq2` | `false` | Skip DESeq2 and all downstream steps |

### Gene ID → Name Conversion

Gene names are extracted automatically from the GTF (`gene_name` attribute). The parameters below are only needed to override with a custom mapping file:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotation` | `null` | Tab-delimited gene ID → name mapping file (gene IDs in col 1). Overrides GTF auto-extraction. |
| `--annotation_key` | `"Gene name"` | Column name in the mapping file that contains gene names |

### Heatmap

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_pheatmap` | `false` | Skip the clustered heatmap step |

---

## Samplesheet Format

A CSV with at minimum `sample` and condition columns. **List control samples first** — the first condition in the samplesheet is used as the DESeq2 reference/control group by default:

```csv
sample,condition,replicate
WT_rep1,WT,1
WT_rep2,WT,2
WT_rep3,WT,3
DKO_rep1,DKO,1
DKO_rep2,DKO,2
DKO_rep3,DKO,3
```

The `sample` column must match the FASTQ filename basenames (before `_R1`/`_R2`).

To explicitly set the control group regardless of samplesheet order, use `--control_group WT`.

---

## Control Group and DESeq2 Comparison Direction

`RunDESeq2.R` uses the **first group label** as the control/reference. The pipeline reorders the count matrix columns so that control-group samples come first, ensuring the output shows `log2(experimental/control)` (e.g. `log2(DKO/WT)`).

The control group is determined by:
1. `--control_group WT` on the command line (takes priority), or
2. The first condition listed in the samplesheet (automatic default)

---

## Output

All results are written to `--outdir` (default: `results/`):

```
results/
├── fastqc/
│   ├── raw/                          FastQC reports for raw reads
│   └── trimmed/                      FastQC reports for trimmed reads
├── trimgalore/                       Trim Galore logs
├── star/
│   ├── <sample>/bam/                 Sorted BAM + index
│   └── <sample>/logs/                STAR alignment logs
├── infer_experiment/                 RSeQC strandedness reports (if --bed set)
├── featurecounts/
│   └── *.featureCounts.txt           Per-sample count files
├── deseq2/
│   ├── merged_counts.txt             Combined count matrix (all samples)
│   ├── merged_counts_named.txt       Gene-name matrix (gene names from GTF, or --annotation override)
│   ├── DESeq2_output.txt             Full DESeq2 results + normalised counts
│   ├── DESeq2_output.PCA.pdf         PCA + sd-mean + boxplot
│   ├── DESeq2_output.dispersion.pdf  Dispersion estimates
│   ├── DESeq2_output.pvalue_histogram.pdf
│   ├── DESeq2_output.MA_plot.no_lfcShrink.pdf
│   ├── MA_plot.*.pdf                 MA plot
│   ├── Volcano_plot.*.pdf            Volcano plot
│   ├── sessionInfo.*.txt             R session info
│   └── heatmap/
│       ├── heatmap_input.txt         rlog matrix for significant genes
│       └── DEsig_heatmap.pdf         Clustered heatmap
├── multiqc/
│   └── rnaseq_multiqc_report.html    Aggregate QC report
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

### DESeq2 output columns (`DESeq2_output.txt`)

| Column | Description |
|--------|-------------|
| `Feature.ID` | Gene name (from GTF auto-extraction, or `--annotation` override; Ensembl ID if GTF has no `gene_name`) |
| `baseMean` | Mean normalised count across all samples |
| `log2(exp/ctrl)` | Shrunk log2 fold change (lfcShrink, normal method) |
| `lfcSE` | Standard error of log2FC |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (BH/FDR) |
| `*.norm` columns | Per-sample normalised counts |
| `*.rlog` columns | Per-sample rlog-transformed counts |

---

## Annotation File Format

Gene names are extracted automatically from the GTF. If you want to override with a custom mapping file (e.g. from Ensembl BioMart), provide a tab-delimited file with Ensembl gene IDs in column 1 and gene names in a named column:

```
Gene stable ID	Gene name	Gene type	Gene description
ENSMUSG00000064372	mt-Tp	Mt_tRNA	mitochondrially encoded tRNA proline
ENSMUSG00000064371	mt-Tt	Mt_tRNA	mitochondrially encoded tRNA threonine
ENSMUSG00000064370	mt-Cytb	protein_coding	mitochondrially encoded cytochrome b
```

This file can be exported from Ensembl BioMart. Use `--annotation_key "Gene name"` to specify which column contains the gene symbols.

---

## Building a STAR Index

If you don't have a pre-built STAR index:

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/annotation.gtf \
     --runThreadN 8
```

The STAR version used for index building must match the pipeline's STAR version (2.7.1a). For genomes < 200 Mb, add `--genomeSAindexNbases 11`.

---

## Repository Structure

```
RNAseq-pipeline/
├── main.nf                     Main Nextflow workflow
├── nextflow.config             All params and profiles
├── environment.yml             Conda environment (includes R + Bioconductor)
├── README.md
├── modules/
│   ├── fastqc.nf
│   ├── trim_galore.nf
│   ├── star_align.nf
│   ├── infer_experiment.nf     RSeQC strandedness auto-detection (all BAMs)
│   ├── featurecounts.nf
│   ├── deseq2.nf               merge + reorder + [annotate] + RunDESeq2.R
│   ├── pheatmap.nf             cluster_draw_pheatmap.R
│   └── multiqc.nf
├── bin/                        Scripts auto-added to $PATH by Nextflow
│   ├── RunDESeq2.R             DESeq2 differential expression
│   ├── draw_MA_plot_from_DESeq2_analysis.R    MA plot
│   ├── draw_volcano_plot_from_DESeq2_analysis.R  Volcano plot
│   ├── cluster_draw_pheatmap.R Heatmap (Cluster 3.0 + ComplexHeatmap)
│   ├── id_to_name.py           Gene ID → name conversion
│   └── merge_counts.py         Merge per-sample featureCounts files
├── conf/
│   ├── base.config             CPU/memory/time per process label
│   ├── slurm.config            SLURM executor settings (partition 24, wibrusers)
│   └── test.config             Lightweight test profile
├── singularity_cache/
│   ├── Singularity.def         Container definition (build with --fakeroot)
│   ├── environment.yml         Environment spec used inside the container
│   └── rnaseqpipe-portable-2.0.sif  (built by user — not in repo)
├── assets/
│   └── multiqc_config.yml
└── test/
    └── samplesheet.csv         Example samplesheet
```

---

## Notes

- MA and Volcano plots are generated by `draw_MA_plot_from_DESeq2_analysis.R` and `draw_volcano_plot_from_DESeq2_analysis.R` in `bin/`. These run in all environments including Singularity.
- Always use **absolute paths** for all file parameters to avoid path resolution issues.
- The `work/` directory (Nextflow cache) is separate from `--outdir`. Do not delete `work/` between runs if you want to use `-resume`.
- To submit the Nextflow head process itself to SLURM (recommended for long runs), wrap the `nextflow run` command in an `srun` call as shown above.

---

## Author

Kalpani De Silva (kdesilva@wi.mit.edu)
Bioinformatics and Research Computing (BaRC), Whitehead Institute for Biomedical Research
Version 2.0-portable, May 2026
