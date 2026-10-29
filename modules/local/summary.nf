process SUMMARY {
    label 'process_single'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path fastp,        stageAs: 'fastp/*'
    path star,         stageAs: 'star/*'
    path assignment,   stageAs: 'in/assignment/*'
    path viral,        stageAs: 'in/viral/*'
    path strandedness, stageAs: 'in/strand/*'
    path deseq,        stageAs: 'deseq/*'
    path fgsea,        stageAs: 'fgsea/*'

    output:
    path 'summary.tsv',   emit: samples
    path 'contrasts.tsv', emit: contrasts

    script:
    """
    collect_summary.py \\
        --fastp fastp/* --star star/* \\
        --assignment ${assignment} --viral ${viral} --strandedness ${strandedness} \\
        --deseq deseq/* --fgsea fgsea/*
    """

    stub:
    """
    touch summary.tsv contrasts.tsv
    """
}
