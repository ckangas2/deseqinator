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
    multiqc --force --filename multiqc_report.html in
    """

    stub:
    """
    touch multiqc_report.html; mkdir multiqc_data
    """
}
