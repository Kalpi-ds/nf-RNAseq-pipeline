#!/usr/bin/env nextflow
////////////////////////////////////////////////////////////////////////////////
//  RNAseq-pipeline — main.nf                                                 //
//  Master workflow: QC → Trim → Align → Quantify → DESeq2 → MultiQC         //
//                                                                            //
//  Usage (paired-end):                                                       //
//    nextflow run main.nf                                                    //
//      --reads "data/*_{R1,R2}.fastq.gz"                                     //
//      --gtf        /path/to/annotation.gtf                                  //
//      --star_index /path/to/star_index                                      //
//      --bed        /path/to/annotation.bed    (for auto strandedness)       //
//      --samplesheet samplesheet.csv                                         //
////////////////////////////////////////////////////////////////////////////////

nextflow.enable.dsl = 2

// -------------------------------------------------------------------------- //
// Import modules
// -------------------------------------------------------------------------- //
include { FASTQC as FASTQC_RAW      } from './modules/fastqc'
include { TRIM_GALORE               } from './modules/trim_galore'
include { FASTQC as FASTQC_TRIMMED  } from './modules/fastqc'
include { STAR_ALIGN                } from './modules/star_align'
include { INFER_EXPERIMENT          } from './modules/infer_experiment'
include { FEATURECOUNTS             } from './modules/featurecounts'
include { DESEQ2                    } from './modules/deseq2'
include { PHEATMAP                  } from './modules/pheatmap'
include { MULTIQC                   } from './modules/multiqc'

// -------------------------------------------------------------------------- //
// Default parameters  (override via CLI or nextflow.config)
// -------------------------------------------------------------------------- //
params.reads             = "*.fastq.gz"
params.single_end        = false
params.outdir            = "results"
params.genome            = null
params.gtf               = null
params.star_index        = null

// Strandedness
// If --bed is provided, strandedness is auto-detected via infer_experiment.py.
// Otherwise, --fc_strandedness is used (0=unstranded, 1=forward, 2=reverse).
params.bed               = null
params.fc_strandedness   = 0

// Trimming
params.min_length        = 20
params.quality_cutoff    = 20
params.skip_trimming     = false

// Alignment
params.threads           = 8
params.star_ram          = 40
params.star_twopass      = false
params.align_intron_min  = 20
params.align_intron_max  = 1000000
params.align_mates_gap   = 1000000

// Quantification
params.fc_group_features     = "gene_id"
params.fc_extra_attributes   = null

// Differential expression
params.samplesheet       = null
params.condition_col     = "condition"
params.control_group     = null
params.lfc_threshold     = 1.0
params.padj_cutoff       = 0.05
params.skip_deseq2       = false

// Gene ID → name conversion (id_to_name.py)
// Provide --annotation pointing to a tab-delimited file with gene IDs in col 1
// and gene names in a column specified by --annotation_key.
// Example: Ensembl gene_id ↔ Gene name mapping from BioMart
params.annotation        = null
params.annotation_key    = "Gene name"

// Heatmap of significant DE genes (cluster_draw_pheatmap.R)
params.skip_pheatmap     = false

// MultiQC
params.multiqc_config    = "${projectDir}/assets/multiqc_config.yml"

// -------------------------------------------------------------------------- //
// Parameter validation
// -------------------------------------------------------------------------- //
if (!params.star_index) error "ERROR: --star_index is required."
if (!params.gtf)        error "ERROR: --gtf is required."
if (!params.skip_deseq2 && !params.samplesheet) {
    log.warn "WARNING: --samplesheet not provided. Skipping DESeq2."
    params.skip_deseq2 = true
}
if (!params.bed) {
    log.warn "WARNING: --bed not provided. Using --fc_strandedness ${params.fc_strandedness}. " +
             "Provide --bed to auto-detect strandedness."
}
if (params.samplesheet) {
    def ss = file(params.samplesheet)
    if (ss.exists()) {
        def lines  = ss.readLines()
        if (lines.size() < 2) error "ERROR: Samplesheet has no data rows: ${params.samplesheet}"
        def header = lines[0].split(',')*.trim()
        if (!header.contains('sample'))
            error "ERROR: Samplesheet missing required 'sample' column: ${params.samplesheet}"
        if (!header.contains(params.condition_col))
            error "ERROR: Samplesheet missing condition column '${params.condition_col}': ${params.samplesheet}"
        def sampleIdx = header.indexOf('sample')
        def samples   = lines[1..-1].collect { it.split(',')[sampleIdx]?.trim() }.findAll { it }
        def dups      = samples.groupBy { it }.findAll { k, v -> v.size() > 1 }.keySet()
        if (dups) error "ERROR: Duplicate sample names in samplesheet: ${dups.join(', ')}"
    }
}

// -------------------------------------------------------------------------- //
// Workflow
// -------------------------------------------------------------------------- //
workflow {

    // -- Input channel -------------------------------------------------------
    if (params.single_end) {
        reads_ch = Channel
            .fromPath(params.reads, checkIfExists: true)
            .map { f -> [ [id: f.simpleName, single_end: true], [ f ] ] }
    } else {
        reads_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .map { id, files -> [ [id: id, single_end: false], files ] }
    }

    star_index = file(params.star_index, checkIfExists: true)
    gtf        = file(params.gtf,        checkIfExists: true)

    // -- QC on raw reads -----------------------------------------------------
    FASTQC_RAW(reads_ch, 'raw')

    // -- Adapter & quality trimming ------------------------------------------
    if (!params.skip_trimming) {
        TRIM_GALORE(reads_ch)
        trimmed_ch        = TRIM_GALORE.out.reads
        trimmed_logs_ch   = TRIM_GALORE.out.log
        FASTQC_TRIMMED(trimmed_ch, 'trimmed')
        trimmed_fastqc_ch = FASTQC_TRIMMED.out.zip
    } else {
        trimmed_ch        = reads_ch
        trimmed_logs_ch   = Channel.empty()
        trimmed_fastqc_ch = Channel.empty()
    }

    // -- Alignment -----------------------------------------------------------
    STAR_ALIGN(trimmed_ch, star_index, gtf)

    // -- Strandedness detection ----------------------------------------------
    // If --bed provided: run infer_experiment.py on ALL BAMs, collect the
    // per-sample strandedness integers, and take the majority vote.
    // Otherwise: fall back to --fc_strandedness parameter.
    if (params.bed) {
        bed_file = file(params.bed, checkIfExists: true)
        INFER_EXPERIMENT(STAR_ALIGN.out.bam, bed_file)
        strandedness_ch = INFER_EXPERIMENT.out.strandedness
            .map  { f -> f.text.trim().toInteger() }
            .collect()
            .map { values ->
                def cnt = [0:0, 1:0, 2:0]
                values.each { v -> cnt[v] = (cnt[v] ?: 0) + 1 }
                def winner = cnt.max { it.value }.key
                log.info "Strandedness consensus across ${values.size()} sample(s): ${cnt} → using ${winner}"
                winner
            }
            .first()    // value channel — broadcasts to all FEATURECOUNTS jobs
    } else {
        strandedness_ch = Channel.value(params.fc_strandedness)
    }

    // -- Quantification ------------------------------------------------------
    // Combine each [meta, bam] with the single strandedness value
    bam_with_strand = STAR_ALIGN.out.bam.combine(strandedness_ch)
    FEATURECOUNTS(bam_with_strand, gtf)

    // -- Differential expression ---------------------------------------------
    if (!params.skip_deseq2) {
        DESEQ2(
            FEATURECOUNTS.out.counts.map { meta, counts -> counts }.collect(),
            file(params.samplesheet, checkIfExists: true),
            params.condition_col,
            params.control_group ?: '',
            gtf
        )

        // -- Heatmap of significant DE genes ---------------------------------
        if (!params.skip_pheatmap) {
            PHEATMAP(DESEQ2.out.results)
        }
    }

    // -- Aggregate QC report -------------------------------------------------
    multiqc_input_ch = Channel.empty()
        .mix( FASTQC_RAW.out.zip.map       { meta, zip -> zip } )
        .mix( trimmed_fastqc_ch.map        { meta, zip -> zip } )
        .mix( trimmed_logs_ch.map          { meta, log -> log } )
        .mix( STAR_ALIGN.out.log_final.map { meta, log -> log } )
        .mix( FEATURECOUNTS.out.summary.map{ meta, s   -> s   } )
        .collect()

    multiqc_config = file(params.multiqc_config, checkIfExists: false)
    MULTIQC(multiqc_input_ch, multiqc_config)

    // -- Completion messages -------------------------------------------------
    STAR_ALIGN.out.bam
        .map { meta, bam -> "✓ BAM: ${bam}" }
        .view()

    FEATURECOUNTS.out.counts
        .map { meta, counts -> "✓ Counts: ${counts}" }
        .view()
}
