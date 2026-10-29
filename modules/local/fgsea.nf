process FGSEA {
    tag "${id}"
    label 'process_low'
    publishDir "${params.outdir}/fgsea/${id}", mode: params.publish_dir_mode

    input:
    tuple val(id), path(deseq)
    path  gmt, stageAs: 'in/gmt/*'
    val   gene_sets

    output:
    path "${id}.fgsea.tsv",         emit: results
    path "${id}.fgsea_summary.tsv", emit: summary
    path "${id}.*",                 emit: all

    script:
    """
    run_fgsea.R \\
        --deseq ${deseq} --id ${id} \\
        --gene_sets '${gene_sets}' --gmt ${gmt} \\
        --min_size ${params.fgsea_min_size} --max_size ${params.fgsea_max_size} \\
        --min_coverage ${params.fgsea_min_coverage} \\
        --max_null_hits ${params.fgsea_max_null_hits} \\
        --fdr ${params.fdr} --seed ${params.seed}
    """

    stub:
    """
    touch ${id}.fgsea.tsv
    printf 'id\\tsets_sig\\n${id}\\t0\\n' > ${id}.fgsea_summary.tsv
    """
}
