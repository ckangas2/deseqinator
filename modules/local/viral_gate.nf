process VIRAL_GATE {
    tag "all samples"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    path loads, stageAs: 'loads/*'
    path samplesheet, stageAs: 'in/*'

    output:
    path 'viral_load.tsv', emit: tsv

    script:
    """
    viral_gate.py loads/* \\
        --samplesheet ${samplesheet} \\
        --mock-conditions '${params.mock_conditions}' \\
        --max-mock ${params.max_mock_viral_frac} \\
        --min-infected ${params.min_infected_viral_frac} \\
        --out viral_load.tsv
    """

    stub:
    """
    viral_gate.py loads/* --samplesheet ${samplesheet} \\
        --mock-conditions '${params.mock_conditions}' \\
        --max-mock ${params.max_mock_viral_frac} \\
        --min-infected ${params.min_infected_viral_frac} --out viral_load.tsv
    """
}
