/*
========================================================================================
    PROCESS: INFER_EXPERIMENT
    Tool   : RSeQC infer_experiment.py
    Purpose: Auto-detect library strandedness from each BAM file.
             Runs on every BAM; main.nf collects all results and takes the
             majority vote, then passes the consensus integer to FEATURECOUNTS.

    Strandedness codes (featureCounts -s):
      0 = unstranded
      1 = forward stranded (e.g. Ligation / ScriptSeq)
      2 = reverse stranded (e.g. dUTP / TruSeq; most common)
========================================================================================
*/

process INFER_EXPERIMENT {
    tag "${meta.id} [strandedness check]"
    label 'process_single'

    publishDir "${params.outdir}/infer_experiment", mode: 'copy'

    input:
    tuple val(meta), path(bam)   // first BAM only
    path  bed                    // gene annotation in BED12 format

    output:
    path "${meta.id}.infer_experiment.txt",  emit: report
    path "${meta.id}.strandedness.txt",      emit: strandedness
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run RSeQC infer_experiment.py
    infer_experiment.py \\
        -i ${bam} \\
        -r ${bed} \\
        -s 200000 \\
        > ${meta.id}.infer_experiment.txt 2>&1

    # Parse output and write strandedness integer to strandedness.txt
    python3 - <<'PYEOF'
import re

with open("${meta.id}.infer_experiment.txt") as f:
    content = f.read()

print("--- infer_experiment.py output ---")
print(content)

paired = "PairEnd" in content

if paired:
    # Reverse stranded:  "1+-,1-+,2++,2--"
    # Forward stranded:  "1++,1--,2+-,2-+"
    m_rev = re.search(r'1\\+\\-,1\\-\\+,2\\+\\+,2\\-\\-[^:]*:\\s*(\\d+\\.\\d+)', content)
    m_fwd = re.search(r'1\\+\\+,1\\-\\-,2\\+\\-,2\\-\\+[^:]*:\\s*(\\d+\\.\\d+)', content)
else:
    # Forward stranded:  "++,--"
    # Reverse stranded:  "+-,-+"
    m_fwd = re.search(r'\\+\\+,\\-\\-[^:]*:\\s*(\\d+\\.\\d+)', content)
    m_rev = re.search(r'\\+\\-,\\-\\+[^:]*:\\s*(\\d+\\.\\d+)', content)

frac_rev = float(m_rev.group(1)) if m_rev else 0.0
frac_fwd = float(m_fwd.group(1)) if m_fwd else 0.0

if frac_rev > 0.6:
    strand = 2
    label  = "reverse stranded (featureCounts -s 2)"
elif frac_fwd > 0.6:
    strand = 1
    label  = "forward stranded (featureCounts -s 1)"
else:
    strand = 0
    label  = "unstranded (featureCounts -s 0)"

print(f"\\nDetected: {label}")
print(f"  frac_reverse = {frac_rev:.3f}")
print(f"  frac_forward = {frac_fwd:.3f}")

with open("${meta.id}.strandedness.txt", "w") as out:
    out.write(str(strand))
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( python3 -c "import RSeQC; print(RSeQC.__version__)" )
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    """
    echo "Stub: assuming reverse stranded (2)" > ${meta.id}.infer_experiment.txt
    echo "2" > ${meta.id}.strandedness.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: 5.0.1
        python: 3.11.0
    END_VERSIONS
    """
}
