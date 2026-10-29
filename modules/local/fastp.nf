process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: params.publish_dir_mode,
        pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads, stageAs: 'in/r?/*')   // R1 and R2 may share a basename
    val   design_ok   // ordering only: no compute until the design is valid

    output:
    tuple val(meta), path("${meta.id}_?.trimmed.fastq.gz"), emit: reads
    path  '*.fastp.json',                        emit: json
    path  '*.fastp.html',                        emit: html

    script:
    // `reads` renders as one path (single-end) or two space-separated paths
    // in samplesheet order; split in bash rather than index a Groovy value
    // whose type depends on how many files arrived.
    def p = meta.id
    if (meta.single_end) {
        """
        fastp --in1 ${reads} --out1 ${p}_1.trimmed.fastq.gz \\
              --thread ${task.cpus} --json ${p}.fastp.json --html ${p}.fastp.html \\
              ${params.fastp_args}
        """
    } else {
        """
        set -- ${reads}
        [ \$# -eq 2 ] || { echo "expected 2 FASTQs for ${p}, got \$#" >&2; exit 1; }
        fastp --in1 \$1 --in2 \$2 \\
              --out1 ${p}_1.trimmed.fastq.gz --out2 ${p}_2.trimmed.fastq.gz \\
              --detect_adapter_for_pe \\
              --thread ${task.cpus} --json ${p}.fastp.json --html ${p}.fastp.html \\
              ${params.fastp_args}
        """
    }

    stub:
    def p = meta.id
    def outs = meta.single_end ? "${p}_1.trimmed.fastq.gz" : "${p}_1.trimmed.fastq.gz ${p}_2.trimmed.fastq.gz"
    """
    for f in ${outs}; do echo | gzip > \$f; done
    printf '{"summary":{"before_filtering":{"total_reads":100},"after_filtering":{"total_reads":95}}}' > ${p}.fastp.json
    touch ${p}.fastp.html
    """
}
