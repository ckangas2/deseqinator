#!/usr/bin/env nextflow
/*
 * deseqinator — bulk RNA-seq of virus infection, FASTQ to gene sets.
 *
 *   design:     samplesheet + contrasts validated before any compute
 *   reference:  validate -> host + virus STAR index
 *   per sample: fastp -> STAR (BAM + GeneCounts) -> viral load
 *   run:        strandedness call -> featureCounts -> viral gate
 *   contrast:   DESeq2 (+ label swaps) -> fgsea (+ rank-shuffle control)
 *
 * Every gate stops the run. Each guards a failure that otherwise produces
 * a complete set of outputs with the wrong numbers in them.
 */
nextflow.enable.dsl = 2

include { CHECK_DESIGN       } from './modules/local/check_design'
include { VALIDATE_REFERENCE } from './modules/local/validate_reference'
include { STAR_INDEX         } from './modules/local/star_index'
include { GENE_TABLE         } from './modules/local/gene_table'
include { FASTP              } from './modules/local/fastp'
include { STAR_ALIGN         } from './modules/local/star_align'
include { VIRAL_LOAD         } from './modules/local/viral_load'
include { VIRAL_GATE         } from './modules/local/viral_gate'
include { INFER_STRANDEDNESS } from './modules/local/infer_strandedness'
include { FEATURECOUNTS      } from './modules/local/featurecounts'
include { DESEQ2             } from './modules/local/deseq2'
include { FGSEA              } from './modules/local/fgsea'
include { RECONCILE          } from './modules/local/reconcile'
include { MULTIQC            } from './modules/local/multiqc'
include { SUMMARY            } from './modules/local/summary'

def helpMessage() {
    log.info """
    deseqinator v${workflow.manifest.version}

    Usage:
      nextflow run ckangas2/deseqinator -profile test,docker
      nextflow run ckangas2/deseqinator -profile docker \\
          --input samplesheet.csv --contrasts contrasts.csv \\
          --host_fasta genome.fa --host_gtf genes.gtf --viral_fasta virus.fa

    Input
      --input              CSV: sample,fastq_1,fastq_2,cell_line,condition,replicate[,batch][,deposited_id]
                           fastq_2 empty for single-end; relative paths resolve
                           against the samplesheet's folder
      --contrasts          CSV: id,cell_line,reference,target  (cell_line '*' = all lines,
                           with cell_line added to the design)

    Reference
      --host_fasta, --host_gtf, --viral_fasta
      --sjdb_overhang      read length - 1 (default: ${params.sjdb_overhang})
      --star_index + --viral_contigs   reuse an index built earlier (both, same build)

    Gates
      --strandedness       auto | unstranded | forward | reverse (default: ${params.strandedness})
      --min_assigned_frac  featureCounts assignment floor (default: ${params.min_assigned_frac})
      --mock_conditions    conditions that must carry no virus (default: ${params.mock_conditions})
      --max_mock_viral_frac / --min_infected_viral_frac
      --min_replicates     per contrast arm (default: ${params.min_replicates})

    DE / GSEA
      --fdr, --lfc_threshold, --min_count, --label_swaps
      --gene_sets          msigdb:<collection> or a .gmt file (default: ${params.gene_sets})
      --fgsea_min_size, --fgsea_max_size, --fgsea_min_coverage, --fgsea_max_null_hits

    Optional
      --deposited_counts   authors' count matrix; needs a deposited_id column in --input

    Output
      --outdir (default: ${params.outdir})   --save_bam   --save_reference   --skip_multiqc
    """.stripIndent()
}

def isSet(v) {
    return v != null && !(v instanceof Boolean) && !(v.toString().trim() in ['', 'null', 'false', 'true'])
}

def validateParams() {
    def problems = []
    if (!isSet(params.input))     problems << "--input is required"
    if (!isSet(params.contrasts)) problems << "--contrasts is required"
    if (!isSet(params.host_gtf))  problems << "--host_gtf is required (featureCounts counts against it)"
    if (isSet(params.star_index)) {
        if (!isSet(params.viral_contigs))
            problems << "--star_index needs --viral_contigs from the same build"
    } else {
        ['host_fasta', 'viral_fasta'].each { p ->
            if (!isSet(params[p])) problems << "--${p} is required when --star_index is not given"
        }
    }
    if (!(params.strandedness in ['auto', 'unstranded', 'forward', 'reverse']))
        problems << "--strandedness must be auto, unstranded, forward or reverse"
    if (!params.gene_sets.toString().startsWith('msigdb:') && !file(params.gene_sets.toString()).exists())
        problems << "--gene_sets must be msigdb:<collection> or an existing .gmt file"
    if (problems) {
        problems.each { msg -> log.error msg }
        error "Fix the parameters above (see --help)."
    }
}

def resolvePath(base, String p) {
    return (p.contains('://') || p.startsWith('/')) ? p : "${base}/${p}"
}

// ---------------------------------------------------------------- workflow
workflow {
    if (params.help) {
        helpMessage()
    } else {
        validateParams()
        main_flow()
    }
}

workflow main_flow {
    main:
    def samplesheet = file(params.input, checkIfExists: true)
    def contrasts   = file(params.contrasts, checkIfExists: true)
    def host_gtf    = file(params.host_gtf, checkIfExists: true)
    def base        = samplesheet.parent

    def rows = samplesheet.splitCsv(header: true, strip: true)
    if (rows.size() == 0) error "samplesheet ${params.input} has a header but no samples"
    def layouts = rows.collect { row -> (row.fastq_2 ?: '') as boolean }.unique()
    if (layouts.size() > 1) error "samplesheet mixes single-end and paired-end libraries; split the run"
    def single_end = !layouts[0]

    // Design first: a bad contrast is caught in seconds, not after alignment.
    CHECK_DESIGN(samplesheet, contrasts)

    def reads = channel.fromList(rows).map { row ->
        def meta = [id: row.sample, single_end: single_end, condition: row.condition,
                    cell_line: row.cell_line]
        def fq = [file(resolvePath(base, row.fastq_1), checkIfExists: true)]
        if (!single_end) fq << file(resolvePath(base, row.fastq_2), checkIfExists: true)
        tuple(meta, fq)
    }

    // ------------------------------------------------ reference
    def index         = null
    def viral_contigs = null
    if (isSet(params.star_index)) {
        index         = channel.value(file(params.star_index, checkIfExists: true))
        viral_contigs = channel.value(file(params.viral_contigs, checkIfExists: true))
    } else {
        def host_fasta  = file(params.host_fasta,  checkIfExists: true)
        def viral_fasta = file(params.viral_fasta, checkIfExists: true)
        VALIDATE_REFERENCE(host_fasta, host_gtf, viral_fasta)
        STAR_INDEX(host_fasta, host_gtf, viral_fasta, VALIDATE_REFERENCE.out.ok)
        index         = STAR_INDEX.out.index
        viral_contigs = STAR_INDEX.out.viral_contigs
    }
    GENE_TABLE(host_gtf)

    // ------------------------------------------------ per sample
    FASTP(reads, CHECK_DESIGN.out.contrasts.map { _f -> true })
    STAR_ALIGN(FASTP.out.reads, index)
    VIRAL_LOAD(STAR_ALIGN.out.bam, viral_contigs)

    // ------------------------------------------------ run-level gates
    VIRAL_GATE(VIRAL_LOAD.out.tsv.collect(), samplesheet)
    INFER_STRANDEDNESS(STAR_ALIGN.out.reads_per_gene.map { _meta, f -> f }.collect())
    def strand = INFER_STRANDEDNESS.out.call.map { f -> f.text.trim() }

    FEATURECOUNTS(
        STAR_ALIGN.out.bam.map { _meta, bam, _bai -> bam }.collect(),
        host_gtf, strand, single_end)

    // ------------------------------------------------ per contrast
    def contrast_ids = CHECK_DESIGN.out.contrasts
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.id }

    DESEQ2(contrast_ids, FEATURECOUNTS.out.matrix, samplesheet,
           CHECK_DESIGN.out.contrasts, GENE_TABLE.out.tsv, VIRAL_GATE.out.tsv)

    def from_msigdb = params.gene_sets.toString().startsWith('msigdb:')
    def gmt = from_msigdb ? file("${projectDir}/assets/NO_FILE") : file(params.gene_sets, checkIfExists: true)
    FGSEA(DESEQ2.out.results, gmt, from_msigdb ? params.gene_sets : 'gmt')

    // ------------------------------------------------ optional
    if (isSet(params.deposited_counts)) {
        RECONCILE(FEATURECOUNTS.out.matrix, GENE_TABLE.out.tsv,
                  file(params.deposited_counts, checkIfExists: true), samplesheet)
    }
    if (!params.skip_multiqc) {
        MULTIQC(FASTP.out.json
                    .mix(STAR_ALIGN.out.log, FEATURECOUNTS.out.summary)
                    .collect())
    }

    SUMMARY(FASTP.out.json.collect(), STAR_ALIGN.out.log.collect(),
            FEATURECOUNTS.out.qc, VIRAL_GATE.out.tsv, INFER_STRANDEDNESS.out.tsv,
            DESEQ2.out.summary.collect(), FGSEA.out.summary.collect())
    SUMMARY.out.contrasts.subscribe { _f ->
        log.info "deseqinator: per-sample QC in ${params.outdir}/summary.tsv, per-contrast results in ${params.outdir}/contrasts.tsv"
    }
}
