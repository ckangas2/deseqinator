process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_star_align'
    publishDir "${params.outdir}/star", mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('.bam') || f.endsWith('.bai') ? (params.save_bam ? f : null) : f }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"),             emit: reads_per_gene
    path  "${meta.id}.Log.final.out",                                     emit: log

    script:
    def p = meta.id
    def sort_ram = task.memory.toBytes().intdiv(2)
    // Mate order fixed by name, not by staging order.
    def fq = meta.single_end ? "${p}_1.trimmed.fastq.gz" : "${p}_1.trimmed.fastq.gz ${p}_2.trimmed.fastq.gz"
    """
    set -euo pipefail
    # Coordinate sorting opens many temporary files per thread.
    ulimit -n 4096 2>/dev/null || true

    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${index} \\
         --readFilesIn ${fq} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${p}. \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMattributes NH HI AS nM \\
         --outSAMunmapped None \\
         --quantMode GeneCounts \\
         --limitBAMsortRAM ${sort_ram}

    mv ${p}.Aligned.sortedByCoord.out.bam ${p}.bam
    samtools index -@ ${task.cpus} ${p}.bam
    """

    stub:
    def p = meta.id
    """
    touch ${p}.bam ${p}.bam.bai
    printf 'N_unmapped\\t0\\t0\\t0\\nN_multimapping\\t0\\t0\\t0\\nN_noFeature\\t0\\t0\\t0\\nN_ambiguous\\t0\\t0\\t0\\nG1\\t1000\\t50\\t950\\n' > ${p}.ReadsPerGene.out.tab
    printf '   Uniquely mapped reads %% |\\t90.00%%\\n' > ${p}.Log.final.out
    """
}
