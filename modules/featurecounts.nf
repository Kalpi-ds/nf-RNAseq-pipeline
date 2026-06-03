/*
========================================================================================
    PROCESS: FEATURECOUNTS
    Tool   : featureCounts (Subread v2.0.3+)
    Purpose: Gene-level read quantification from sorted BAM files.
             Strandedness (-s) is passed as a channel value, auto-detected by
             INFER_EXPERIMENT or set manually via --fc_strandedness.
========================================================================================
*/

process FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/featurecounts/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), val(strandedness)   // strandedness: 0, 1, or 2
    path  gtf                                        // GTF annotation

    output:
    tuple val(meta), path("*.featureCounts.txt"),           emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"),   emit: summary
    path  "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def extra_args     = task.ext.args ?: ''
    def group_features = params.fc_group_features
    def extra_attr     = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
    def pe_flags       = meta.single_end ? '' : '-p'

    """
    featureCounts \\
        -a ${gtf} \\
        -o ${prefix}.featureCounts.txt \\
        -g ${group_features} \\
        -T ${task.cpus} \\
        -s ${strandedness} \\
        ${pe_flags} \\
        ${extra_attr} \\
        ${extra_args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( featureCounts -v 2>&1 | grep -oP '(?<=version ).*' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.txt ${prefix}.featureCounts.txt.summary
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: 2.0.3
    END_VERSIONS
    """
}
