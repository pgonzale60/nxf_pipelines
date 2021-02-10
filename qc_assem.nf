nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "multiassem-${date}"
params.reads = "mini_DF5120.ccs.fasta.gz"
params.assemblies = "mini_DF5120.hifiasm.fasta.gz"
params.dmnd_db = "nigons.dmnd"
params.max_target_seqs = 1000
params.evalue = 0.0001


reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.telo)?(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }



assemblies = Channel
                .fromPath(params.assemblies, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.hifiasm)?(\.flye)?(\.wtdbg2)?(\.canu)?(\.fa)?(\.fasta)?(\.gz)?$/, file.Name - ~/(\.fa)?(\.fasta)?(\.gz)?$/, file) }

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()


process mask_assembly {
    tag "${assembler}"
    label 'btk'

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(strain), val(assembler), path("${assembler}.masked.fasta")

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $assembly > assembly.fasta
        else
            ln -s $assembly assembly.fasta
      fi
      windowmasker -in assembly.fasta \
                      -infmt fasta \
                      -mk_counts \
                      -sformat obinary \
                      -out tmp.counts 2> log \
        && windowmasker -in assembly.fasta \
                        -infmt fasta \
                        -ustat tmp.counts \
                        -dust T \
                        -outfmt fasta \
                        -out ${assembler}.masked.fasta
     rm assembly.fasta
      """
}


process chunk_assembly {
    tag "${assembler}"
    label 'btk'

    input:
      tuple val(strain), val(assembler), path(assembly)

    output:
      tuple val(strain), val(assembler), path("${assembler}.chunks.fasta")

    script:
      """
      chunk_fasta.py --in ${assembly} \
        --chunk 100000 --overlap 0 --max-chunks 20 \
        --out ${assembler}.chunks.fasta
      """
}

process diamond_search {
    tag "${assembler}"
    label 'btk'

    input:
      tuple val(strain), val(assembler), path(assembly)
      path(dmnd_db)

    output:
      tuple val(strain), val(assembler), path("${assembler}.chunk.diamond.tsv")

    script:
      """
      diamond blastx \
            --query ${assembly} \
            --db $params.dmnd_db \
            --outfmt 6 qseqid qstart qend sseqid sstart send bitscore scovhsp \
            --max-target-seqs $params.max_target_seqs \
            --max-hsps 1 \
            --evalue $params.evalue \
            --threads ${task.cpus} \
            > ${assembler}.chunk.diamond.tsv
      """
}

process unchunk_hits {
    tag "${assembler}"
    publishDir "$params.outdir/nigons", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), val(assembler), path(diamondHits)

    output:
      tuple val(strain), val(assembler), path("${assembler}.diamond.tsv.gz")

    script:
      """
      unchunk_ng_dmnd.py --in $diamondHits \
        --out ${assembler}.diamond.tsv
        gzip ${assembler}.diamond.tsv
      """
}

process map_telomeric_reads {
    tag "${assembler}"
    publishDir "$params.outdir/teloMaps", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads), val(assembler), path(assembly)

    output:
      tuple val(strain), val(assembler), path("${assembler}.teloMapped.mmap2.tsv.gz")

    script:
      """
      minimap2 ${assembly} $reads | \
        gzip -c > ${assembler}.teloMapped.mmap2.tsv.gz
      """
}



workflow {
    mask_assembly(assemblies) | chunk_assembly
    diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
    map_telomeric_reads(reads.join(assemblies))
}
