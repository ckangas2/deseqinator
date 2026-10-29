process GENE_TABLE {
    tag "${gtf.name}"
    label 'process_single'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path gtf, stageAs: 'in/*'

    output:
    path 'gene_annotation.tsv', emit: tsv

    script:
    """
    gtf_gene_table.py --gtf ${gtf} --out gene_annotation.tsv
    """

    stub:
    """
    gtf_gene_table.py --gtf ${gtf} --out gene_annotation.tsv
    """
}
