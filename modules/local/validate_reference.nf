process VALIDATE_REFERENCE {
    tag "host+virus"
    label 'process_single'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path host_fasta,  stageAs: 'in/host/*'
    path host_gtf,    stageAs: 'in/host_gtf/*'
    path viral_fasta, stageAs: 'in/viral/*'

    output:
    path 'validate_reference.log', emit: log
    val  true,                     emit: ok

    script:
    """
    set -euo pipefail
    zcat -f ${host_fasta} | grep '^>' > host_contigs.txt || true
    validate_reference.py \\
        --host-contigs host_contigs.txt \\
        --host-gtf ${host_gtf} \\
        --viral-fasta ${viral_fasta} \\
        > validate_reference.log 2>&1 || { cat validate_reference.log >&2; exit 1; }
    cat validate_reference.log
    """

    stub:
    """
    echo "stub" > validate_reference.log
    """
}
