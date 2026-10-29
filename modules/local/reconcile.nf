process RECONCILE {
    tag "vs deposited"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    path counts,      stageAs: 'in/counts/*'
    path annotation,  stageAs: 'in/annotation/*'
    path deposited,   stageAs: 'in/deposited/*'
    path samplesheet, stageAs: 'in/samplesheet/*'

    output:
    path 'reconciliation.tsv',        emit: tsv
    path 'reconciliation_matrix.tsv', emit: matrix

    script:
    """
    reconcile_counts.py --ours ${counts} --annotation ${annotation} \\
        --deposited ${deposited} --samplesheet ${samplesheet} \\
        --min-rho ${params.reconcile_min_rho}
    """

    stub:
    """
    touch reconciliation.tsv reconciliation_matrix.tsv
    """
}
