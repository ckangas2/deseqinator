process CHECK_DESIGN {
    tag "${samplesheet.name}"
    label 'process_single'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path samplesheet, stageAs: 'in/samplesheet/*'
    path contrasts,   stageAs: 'in/contrasts/*'

    output:
    path 'design/contrasts.validated.tsv', emit: contrasts
    path 'design/design_report.txt',       emit: report

    script:
    """
    check_design.py \\
        --samplesheet ${samplesheet} \\
        --contrasts ${contrasts} \\
        --min-replicates ${params.min_replicates} \\
        --out-dir design
    """

    stub:
    """
    check_design.py --samplesheet ${samplesheet} --contrasts ${contrasts} \\
        --min-replicates ${params.min_replicates} --out-dir design
    """
}
