/*
========================================================================================
    PROCESSES: STAR_GENOMEGENERATE  &  STAR_ALIGN
    Tool     : STAR v2.7.1a
    Purpose  : (1) Build genome index from FASTA + GTF
               (2) Align single-end or paired-end reads to reference genome
========================================================================================
*/

// ---------------------------------------------------------------------------
// Build STAR genome index
// ---------------------------------------------------------------------------
process STAR_GENOMEGENERATE {
    tag "Building STAR index"
    label 'process_high'

    publishDir "${params.outdir}/genome/star_index", mode: 'copy'

    input:
    path fasta   // Reference genome FASTA (e.g. GRCh38.primary_assembly.genome.fa)
    path gtf     // Gene annotation GTF

    output:
    path "star_index",   emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genomeSAindexNbases = 14   // Default; reduce to 11 for small genomes
    def extra_args = task.ext.args ?: ''
    """
    mkdir -p star_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases ${genomeSAindexNbases} \\
        --runThreadN ${task.cpus} \\
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$( STAR --version | sed -e "s/STAR_//g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index && touch star_index/Genome
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.1a
    END_VERSIONS
    """
}

// ---------------------------------------------------------------------------
// Align reads with STAR
// ---------------------------------------------------------------------------
process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/star/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.bam'))    "bam/${filename}"
            else if (filename.endsWith('.bai')) "bam/${filename}"
            else if (filename.endsWith('Log.final.out')) "logs/${filename}"
            else if (filename.endsWith('Log.out'))       "logs/${filename}"
            else if (filename.endsWith('SJ.out.tab'))    "logs/${filename}"
            else null
        }

    input:
    tuple val(meta), path(reads)   // [meta, reads] — single-end: [R1.fq.gz], paired-end: [R1.fq.gz, R2.fq.gz]
    path  index                    // STAR index directory
    path  gtf                      // GTF annotation

    output:
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam"),         emit: bam
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam.bai"),     emit: bai,     optional: true
    tuple val(meta), path("*Log.final.out"),                          emit: log_final
    tuple val(meta), path("*Log.out"),                                emit: log_out
    tuple val(meta), path("*SJ.out.tab"),                             emit: sj
    tuple val(meta), path("*ReadsPerGene.out.tab"),                   emit: read_per_gene, optional: true
    path  "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def extra_args   = task.ext.args ?: ''
    def twopass      = params.star_twopass ? '--twopassMode Basic' : ''
    def reads_input  = meta.single_end ? "--readFilesIn ${reads[0]}" : "--readFilesIn ${reads[0]} ${reads[1]}"
    def mates_gap    = meta.single_end ? '' : "--alignMatesGapMax ${params.align_mates_gap}"
    def out_filter_mismatch = 999   // Allow any mismatches (featureCounts handles NH filter)
    """
    # Run STAR alignment
    STAR \\
        --runMode alignReads \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        ${reads_input} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --alignIntronMin ${params.align_intron_min} \\
        --alignIntronMax ${params.align_intron_max} \\
        ${mates_gap} \\
        --outFilterMismatchNmax ${out_filter_mismatch} \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --quantMode GeneCounts \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        --limitBAMsortRAM ${task.memory.toBytes()} \\
        ${twopass} \\
        ${extra_args}

    # Index BAM with samtools
    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$( STAR --version | sed -e "s/STAR_//g" )
        samtools: \$( samtools --version | head -1 | sed 's/samtools //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam.bai
    touch ${prefix}.Log.final.out ${prefix}.Log.out ${prefix}.SJ.out.tab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.1a
        samtools: 1.18
    END_VERSIONS
    """
}
