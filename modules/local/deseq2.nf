process DESEQ2 {
    tag "${id}"
    label 'process_deseq2'
    // Outputs go into a folder named after the contrast; publishDir names only
    // the parent, since Nextflow 26 resolves directives without task inputs.
    publishDir "${params.outdir}/deseq2", mode: params.publish_dir_mode

    input:
    val  id
    path counts,      stageAs: 'in/counts/*'
    path samplesheet, stageAs: 'in/samplesheet/*'
    path contrasts,   stageAs: 'in/contrasts/*'
    path annotation,  stageAs: 'in/annotation/*'
    path viral_gate   // ordering only: no DE on samples that failed the gate

    output:
    tuple val(id), path("${id}/${id}.deseq2.tsv"), emit: results
    path "${id}/${id}.deseq2_summary.tsv",         emit: summary
    path "${id}/*",                                emit: all

    script:
    """
    mkdir ${id} && cd ${id}
    run_deseq2.R \\
        --counts ../${counts} --samplesheet ../${samplesheet} --contrasts ../${contrasts} \\
        --id ${id} --annotation ../${annotation} \\
        --fdr ${params.fdr} --lfc ${params.lfc_threshold} --min_count ${params.min_count} \\
        --label_swaps ${params.label_swaps} --seed ${params.seed}
    """

    stub:
    """
    mkdir ${id}
    printf 'gene_id\\tgene_name\\tstat\\n' > ${id}/${id}.deseq2.tsv
    printf 'id\\tn_up\\tn_down\\n${id}\\t0\\t0\\n' > ${id}/${id}.deseq2_summary.tsv
    """
}
