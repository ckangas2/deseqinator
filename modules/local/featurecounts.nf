process FEATURECOUNTS {
    tag "${bams.size()} samples, ${strandedness}"
    label 'process_medium'
    publishDir "${params.outdir}/counts", mode: params.publish_dir_mode

    input:
    path bams, stageAs: 'bam/*'
    path gtf,  stageAs: 'in/*'
    val  strandedness
    val  single_end

    output:
    path 'counts.tsv',                    emit: matrix
    path 'assignment.tsv',                emit: qc
    path 'featurecounts.txt.summary',     emit: summary
    path 'featurecounts.txt',             emit: raw

    script:
    def s = [unstranded: 0, forward: 1, reverse: 2][strandedness]
    if (s == null) error "unknown strandedness '${strandedness}'"
    def pe = single_end ? '' : '-p --countReadPairs'
    """
    set -euo pipefail
    zcat -f ${gtf} > host.gtf
    # Run from bam/ so the column names are <sample>.bam.
    cd bam
    featureCounts -T ${task.cpus} -a ../host.gtf -o ../featurecounts.txt \\
        -t exon -g gene_id -s ${s} ${pe} *.bam
    cd ..
    rm host.gtf
    tidy_counts.py --featurecounts featurecounts.txt \\
        --min-assigned ${params.min_assigned_frac} \\
        --out-matrix counts.tsv --out-qc assignment.tsv
    """

    stub:
    """
    echo stub > featurecounts.txt
    echo stub > featurecounts.txt.summary
    printf 'sample\\ttotal\\tassigned\\tassigned_frac\\ttop_unassigned\\n' > assignment.tsv
    touch counts.tsv
    """
}
