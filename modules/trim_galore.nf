/*
========================================================================================
    PROCESS: TRIM_GALORE
    Tool   : Trim Galore v0.6.10+ (wraps Cutadapt + FastQC)
    Purpose: Adapter trimming and quality filtering of single-end or paired-end reads
========================================================================================
*/

process TRIM_GALORE {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/trimgalore/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.html')) "fastqc/${filename}"
            else if (filename.endsWith('.zip'))  "fastqc/${filename}"
            else if (filename.endsWith('trimming_report.txt')) "logs/${filename}"
            else null  // don't publish trimmed FASTQs to save disk; remove null to keep them
        }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{_1,_2}.fq.gz"),         emit: reads
    tuple val(meta), path("*trimming_report.txt"),    emit: log
    tuple val(meta), path("*.html"),                  emit: html,     optional: true
    tuple val(meta), path("*.zip"),                   emit: zip,      optional: true
    path  "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cores        = task.cpus > 4 ? 4 : task.cpus    // Trim Galore caps at 4 cores
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def extra_args   = task.ext.args ?: ''
    def paired_flag  = meta.single_end ? '' : '--paired'
    def input_reads  = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    def rename_cmd   = meta.single_end
        ? "mv ${prefix}_trimmed.fq.gz ${prefix}_1.fq.gz"
        : "mv ${prefix}_val_1.fq.gz ${prefix}_1.fq.gz && mv ${prefix}_val_2.fq.gz ${prefix}_2.fq.gz"
    """
    trim_galore \\
        --cores ${cores} \\
        ${paired_flag} \\
        --quality 20 \\
        --length 20 \\
        --fastqc \\
        --gzip \\
        --basename ${prefix} \\
        ${extra_args} \\
        ${input_reads}

    # Rename output to predictable _1 (and _2 for paired-end) pattern
    ${rename_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$( trim_galore --version | grep 'version' | sed 's/.*version //; s/ .*//' )
        cutadapt: \$( cutadapt --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1.fq.gz
    if [ "${meta.single_end}" != "true" ]; then touch ${prefix}_2.fq.gz; fi
    touch ${prefix}_R1_trimming_report.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: 0.6.10
        cutadapt: 4.4
    END_VERSIONS
    """
}
