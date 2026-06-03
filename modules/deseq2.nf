/*
========================================================================================
    PROCESS: DESEQ2
    Tool   : R / DESeq2 (Bioconductor) — BaRC RunDESeq2.R
    Purpose: Merge per-sample featureCounts files, optionally convert gene IDs to
             gene names (id_to_name.py), then run differential expression analysis.

    Script interfaces
    -----------------
    merge_counts.py   --input *.featureCounts.txt --output merged_counts.txt
    id_to_name.py     -r <ref.txt> -q merged_counts.txt -k 'Gene name' > named.txt
    RunDESeq2.R       inputCounts outputFile cond_samp1 cond_samp2 ...
                        (groups are per-sample, in column order of the count matrix)

    Notes
    -----
    • RunDESeq2.R also calls draw_MA_plot/draw_volcano_plot scripts from
      /nfs/BaRC_Public/BaRC_code/R/DESeq2/ — these are available on the cluster NFS
      but not inside the container. The system() calls fail silently; all other output
      (PCA, dispersion, MA-no-shrink PDFs) is still produced.
    • Supply --annotation and --annotation_key to enable gene-name conversion.
========================================================================================
*/

process DESEQ2 {
    label 'process_medium'

    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts_files   // all per-sample .featureCounts.txt files (collected)
    path samplesheet    // CSV: sample,<condition_col>[,batch,...]
    val  condition_col  // column in samplesheet that holds condition labels
    val  control_group  // reference/control condition name (e.g. "WT"); empty = auto
    path gtf            // GTF annotation — used to extract gene_id → gene_name mapping

    output:
    path "merged_counts.txt",         emit: merged_counts
    path "DESeq2_output.txt",         emit: results
    path "merged_counts_named.txt",   emit: named_counts,  optional: true
    path "*.pdf",                     emit: plots,         optional: true
    path "sessionInfo.*.txt",         emit: session_info,  optional: true
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── Step 1: merge per-sample featureCounts files ──────────────────────────
    merge_counts.py \\
        --input *.featureCounts.txt \\
        --output merged_counts.txt

    # ── Step 2: reorder columns so control-group samples come first ───────────
    # RunDESeq2.R uses groups[1] as the control/reference, so column order
    # in the counts matrix determines the comparison direction.
    SAMPLE_GROUPS=\$(python3 - << 'PYEOF'
import csv, sys, pandas as pd

samplesheet_file = "${samplesheet}"
cond_col         = "${condition_col}"
control_group    = "${control_group}"   # empty string if not specified

# Read samplesheet: sample -> condition mapping
sample_cond = {}
with open(samplesheet_file) as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        sample_cond[row["sample"].strip()] = row[cond_col].strip()

# Read the merged counts matrix
df = pd.read_csv("merged_counts.txt", sep="\\t")
samples = list(df.columns[1:])   # drop 'gene_id' / 'Gene.ID'

# Validate: every column must be in the samplesheet
for s in samples:
    if s not in sample_cond:
        sys.exit(
            f"ERROR: sample '{s}' in merged_counts.txt not found in "
            f"samplesheet column 'sample'.  Check your samplesheet."
        )

# Determine control group
if not control_group:
    # Default: use the first condition listed in the samplesheet
    with open(samplesheet_file) as fh:
        reader = csv.DictReader(fh)
        control_group = next(reader)[cond_col].strip()
    print(f"Auto-detected control group: {control_group}", file=sys.stderr)

# Reorder columns: control-group samples first, then the rest
ctrl_samples  = [s for s in samples if sample_cond[s] == control_group]
other_samples = [s for s in samples if sample_cond[s] != control_group]
new_order     = ctrl_samples + other_samples

if not ctrl_samples:
    sys.exit(
        f"ERROR: control_group '{control_group}' not found among conditions: "
        f"{sorted(set(sample_cond[s] for s in samples))}"
    )

# Rewrite merged_counts.txt with reordered columns
df = df[[df.columns[0]] + new_order]
df.to_csv("merged_counts.txt", sep="\\t", index=False)

# Build group labels in the new column order
groups = [sample_cond[s] for s in new_order]
print(" ".join(groups))
PYEOF
)

    # ── Step 3: gene ID → gene name conversion ───────────────────────────────
    # Auto-extracts gene_id → gene_name from GTF unless --annotation is provided.
    # Runs AFTER column reordering so control samples are first in output.
    COUNTS=merged_counts.txt
    ANNOT_REF=""
    ANNOT_KEY="Gene name"

    if [ -n "${params.annotation ?: ''}" ]; then
        ANNOT_REF="${params.annotation}"
        ANNOT_KEY="${params.annotation_key ?: 'Gene name'}"
    else
        python3 - << 'PYEOF'
import re, sys
mapping = {}
with open("${gtf}") as f:
    for line in f:
        if line.startswith('#'):
            continue
        m_id   = re.search(r'gene_id "([^"]+)"', line)
        m_name = re.search(r'gene_name "([^"]+)"', line)
        if m_id and m_name:
            mapping[m_id.group(1)] = m_name.group(1)
if not mapping:
    print("WARNING: No gene_name attribute found in GTF — keeping Ensembl IDs", file=sys.stderr)
    sys.exit(0)
with open("gtf_annotation.txt", "w") as out:
    out.write("gene_id\\tGene name\\n")
    for gid, gname in mapping.items():
        out.write(gid + "\\t" + gname + "\\n")
print(f"Extracted {len(mapping)} gene_id -> gene_name mappings from GTF", file=sys.stderr)
PYEOF
        [ -f gtf_annotation.txt ] && ANNOT_REF="gtf_annotation.txt"
    fi

    if [ -n "\$ANNOT_REF" ]; then
        id_to_name.py -r "\$ANNOT_REF" -q merged_counts.txt -k "\$ANNOT_KEY" > merged_counts_named.txt
        COUNTS=merged_counts_named.txt
    fi

    # ── Step 4: run RunDESeq2.R ───────────────────────────────────────────────
    # SAMPLE_GROUPS is a space-separated string of per-sample condition labels,
    # in the same column order as the count matrix.
    # Note: we avoid the variable name GROUPS because it is a bash built-in
    # (holds the current user's Unix group IDs) and gets silently overwritten.
    # Unset LD_LIBRARY_PATH to prevent conda's shared libraries from
    # interfering with R package loading (can cause segfaults).

    # Pass pipeline thresholds to RunDESeq2.R via environment variables
    export LFC_THRESHOLD="${params.lfc_threshold}"
    export PADJ_CUTOFF="${params.padj_cutoff}"

    SAVED_LD_PATH=\${LD_LIBRARY_PATH:-}
    unset LD_LIBRARY_PATH
    RunDESeq2.R \\
        \$COUNTS \\
        DESeq2_output.txt \\
        \$SAMPLE_GROUPS
    export LD_LIBRARY_PATH=\$SAVED_LD_PATH

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( unset LD_LIBRARY_PATH; R --version | head -1 | sed 's/R version //; s/ .*//' )
        bioc-deseq2: \$( unset LD_LIBRARY_PATH; Rscript -e "cat(as.character(packageVersion('DESeq2')))" )
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    """
    touch merged_counts.txt DESeq2_output.txt
    touch DESeq2_output.PCA.pdf DESeq2_output.dispersion.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.1
        bioc-deseq2: 1.42.0
        python: 3.11.0
    END_VERSIONS
    """
}
