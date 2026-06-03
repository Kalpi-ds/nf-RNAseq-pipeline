/*
========================================================================================
    PROCESS: PHEATMAP
    Tool   : Cluster 3.0 + R cluster_draw_pheatmap.R (BaRC)
    Purpose: Generate a clustered heatmap for significantly DE genes.
             Filters the DESeq2 output for significant genes, extracts rlog-
             transformed counts, and passes the matrix to cluster_draw_pheatmap.R.

    Requirements
    ------------
    • Cluster 3.0 binary at /usr/local/bin/cluster (installed via conda cluster3
      package; a symlink is created in the Docker/Singularity image)
    • R packages: pheatmap, ComplexHeatmap (bioconductor-complexheatmap)

    Parameters (from nextflow.config / CLI)
    ----------------------------------------
    params.padj_cutoff    — FDR threshold for significance  (default 0.05)
    params.lfc_threshold  — |log2FC| threshold              (default 1.0)
    params.skip_pheatmap  — set true to skip this process   (default false)
========================================================================================
*/

process PHEATMAP {
    label 'process_single'

    publishDir "${params.outdir}/deseq2/heatmap", mode: 'copy'

    input:
    path deseq2_results   // DESeq2_output.txt from the DESEQ2 process

    output:
    path "DEsig_heatmap.pdf",    emit: heatmap,       optional: true
    path "heatmap_input.txt",    emit: matrix,        optional: true
    path "*.pdf",                emit: extra_pdfs,    optional: true
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def padj_cut = params.padj_cutoff  ?: 0.05
    def lfc_thr  = params.lfc_threshold ?: 1.0
    """
    # ── Step 1: filter DESeq2 output → heatmap input matrix ──────────────────
    # Keeps genes where padj <= threshold AND |log2FC| >= threshold.
    # Extracts the .rlog columns (rlog-transformed counts) for display.
    python3 - << 'PYEOF'
import sys
import pandas as pd

padj_cut = ${padj_cut}
lfc_thr  = ${lfc_thr}

df = pd.read_csv("${deseq2_results}", sep="\\t")

# Identify the log2FC column (second column, name like "log2(DKO/WT)")
lfc_col = df.columns[1]

# Identify rlog-transformed count columns (end with '.rlog')
rlog_cols = [c for c in df.columns if c.endswith(".rlog")]

if not rlog_cols:
    sys.exit("ERROR: No '.rlog' columns found in DESeq2 output. "
             "Ensure RunDESeq2.R completed successfully.")

# Filter: significant and above fold-change threshold
sig = df[
    (df["padj"].notna()) &
    (df["padj"]  <= padj_cut) &
    (df[lfc_col].abs() >= lfc_thr)
].copy()

if len(sig) == 0:
    print(f"WARNING: No significant genes found "
          f"(padj <= {padj_cut}, |log2FC| >= {lfc_thr}).  "
          f"Heatmap will be skipped.")
    sys.exit(0)

print(f"Heatmap: {len(sig)} significant genes")

# Build output matrix: gene names in col 1, rlog values for each sample
out = sig[["Feature.ID"] + rlog_cols].copy()
# Clean up column names: remove '.rlog' suffix for readability
out.columns = ["Gene"] + [c.replace(".rlog", "") for c in rlog_cols]
out.to_csv("heatmap_input.txt", sep="\\t", index=False)
PYEOF

    # ── Step 2: run cluster_draw_pheatmap.R (only if matrix was created) ─────
    # Unset LD_LIBRARY_PATH to prevent conda's shared libraries from
    # interfering with R package loading (can cause segfaults).
    SAVED_LD_PATH=\${LD_LIBRARY_PATH:-}
    unset LD_LIBRARY_PATH
    if [[ -f heatmap_input.txt ]]; then
        cluster_draw_pheatmap.R heatmap_input.txt DEsig_heatmap.pdf
    fi
    export LD_LIBRARY_PATH=\$SAVED_LD_PATH

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( unset LD_LIBRARY_PATH; R --version | head -1 | sed 's/R version //; s/ .*//' )
        r-complexheatmap: \$( unset LD_LIBRARY_PATH; Rscript -e "cat(as.character(packageVersion('ComplexHeatmap')))" )
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    """
    touch heatmap_input.txt DEsig_heatmap.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.1
        r-complexheatmap: 2.18.0
        python: 3.11.0
    END_VERSIONS
    """
}
