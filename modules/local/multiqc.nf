process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path files, stageAs: 'in/?/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data',        emit: data

    script:
    """
    multiqc --force in    # default names: multiqc_report.html, multiqc_data/
    """

    stub:
    """
    touch multiqc_report.html; mkdir multiqc_data
    """
}
