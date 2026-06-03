/*
========================================================================================
    PROCESS: FASTQC
    Tool   : FastQC v0.11.9+
    Purpose: Quality control on raw or trimmed paired-end FASTQ reads
========================================================================================
*/

process FASTQC {
    tag "${meta.id} [${stage}]"
    label 'process_medium'

    publishDir "${params.outdir}/fastqc/${stage}/${meta.id}", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.zip') ? null : filename }

    input:
    tuple val(meta), path(reads)
    val   stage                   // 'raw' or 'trimmed' — used for output subdirectory

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    fastqc \\
        --threads ${threads} \\
        --outdir . \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_fastqc.html ${prefix}_R1_fastqc.zip
    touch ${prefix}_R2_fastqc.html ${prefix}_R2_fastqc.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.11.9
    END_VERSIONS
    """
}
