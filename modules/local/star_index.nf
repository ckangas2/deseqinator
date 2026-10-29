process STAR_INDEX {
    tag "host+virus"
    label 'process_star_index'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { f -> (f == 'star_index' && !params.save_reference) ? null : f }

    // Inputs staged into their own folders so no output can be written
    // through an input symlink (lesson from aligninator: genome.fa / genes.gtf).
    input:
    path host_fasta,  stageAs: 'in/host/*'
    path host_gtf,    stageAs: 'in/host_gtf/*'
    path viral_fasta, stageAs: 'in/viral/*'
    val  validated    // ordering only: never index an unvalidated reference

    output:
    path 'star_index',           emit: index
    path 'viral_contigs.txt',    emit: viral_contigs
    path 'star_index_check.log', emit: log

    script:
    """
    set -euo pipefail
    set -o noclobber

    # Host first: contig order matches the host-only reference.
    zcat -f ${host_fasta} ${viral_fasta} > merged.genome.fa
    zcat -f ${host_gtf} > host.genes.gtf
    zcat -f ${viral_fasta} | grep '^>' | sed 's/^>//; s/[[:space:]].*//' > viral_contigs.txt

    glen=\$(grep -v '^>' merged.genome.fa | tr -d '\\n' | wc -c)
    sa=\$(python3 -c "import math; print(min(14, int(math.log2(\$glen)/2 - 1)))")

    mkdir star_index
    STAR --runMode genomeGenerate \\
         --runThreadN ${task.cpus} \\
         --genomeDir star_index \\
         --genomeFastaFiles merged.genome.fa \\
         --sjdbGTFfile host.genes.gtf \\
         --sjdbOverhang ${params.sjdb_overhang} \\
         --genomeSAindexNbases \$sa \\
         --limitGenomeGenerateRAM ${task.memory.toBytes()} \\
         --outFileNamePrefix star_index/
    rm -f merged.genome.fa host.genes.gtf

    # The viral contig must be in the index, or viral load reads as zero.
    missing=0
    while read -r c; do
      grep -qx "\$c" star_index/chrName.txt || { echo "viral contig \$c not in index" >&2; missing=1; }
    done < viral_contigs.txt
    n_genes=\$(tail -n +2 star_index/geneInfo.tab | wc -l)
    {
      echo "genome length       : \$glen"
      echo "genomeSAindexNbases : \$sa"
      echo "sjdbOverhang        : ${params.sjdb_overhang}"
      echo "annotated genes     : \$n_genes"
      echo "viral contigs       : \$(tr '\\n' ' ' < viral_contigs.txt)"
    } | tee star_index_check.log
    [ "\$missing" -eq 0 ] || exit 1
    [ "\$n_genes" -gt 0 ] || { echo "FATAL: no annotated genes in the index" >&2; exit 1; }
    """

    stub:
    """
    mkdir star_index
    zcat -f ${viral_fasta} | grep '^>' | sed 's/^>//; s/[[:space:]].*//' > viral_contigs.txt
    echo stub > star_index_check.log
    """
}
