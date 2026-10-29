process VIRAL_LOAD {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)
    path  viral_contigs

    output:
    path "${meta.id}.viral.tsv", emit: tsv

    script:
    // Primary alignments only (-F 0x904: not unmapped, secondary or
    // supplementary), so multimappers count once.
    """
    set -euo pipefail
    total=\$(samtools view -c -F 0x904 ${bam})
    viral=\$(samtools view -c -F 0x904 ${bam} \$(tr '\\n' ' ' < ${viral_contigs}))
    printf '%s\\t%s\\t%s\\n' ${meta.id} "\$total" "\$viral" > ${meta.id}.viral.tsv
    """

    stub:
    def viral = meta.condition in params.mock_conditions.tokenize(',') ? 0 : 300
    """
    printf '%s\\t%s\\t%s\\n' ${meta.id} 1000 ${viral} > ${meta.id}.viral.tsv
    """
}
