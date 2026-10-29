process INFER_STRANDEDNESS {
    tag "all samples"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    path tabs, stageAs: 'rpg/*'

    output:
    path 'strandedness.tsv', emit: tsv
    path 'strandedness.txt', emit: call

    script:
    """
    infer_strandedness.py rpg/* \\
        --declared ${params.strandedness} \\
        --threshold ${params.strand_threshold} \\
        --out-tsv strandedness.tsv --out-call strandedness.txt
    """

    stub:
    """
    infer_strandedness.py rpg/* --declared ${params.strandedness} \\
        --threshold ${params.strand_threshold} \\
        --min-reads 100 --out-tsv strandedness.tsv --out-call strandedness.txt
    """
}
