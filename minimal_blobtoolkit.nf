nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "miniBtk-${date}"
params.reads = "/home/ubuntu/oscheius/0-inputs/test.ccsf.fasta.gz"
params.assemblies = "/home/ubuntu/oscheius/0-inputs/test.fasta"
params.btkPath = "~/sw/blobtoolkit/"
params.blobtoolsPath = "${params.btkPath}/blobtools2/blobtools"
params.dmnd_db = "custom.dmnd"
params.max_target_seqs = 1000
params.evalue = 0.000001
params.taxid = 2613844
params.taxdump = "${params.btkPath}/taxdump/"
params.taxrule = "bestsumorder"




// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/067/145/GCA_904067145.1_BOKI2/GCA_904067145.1_BOKI2_protein.faa.gz # Bursaphelenchus okinawaensis
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz # Caenorhabditis elegans
// wget .../GCF_000005845.2_ASM584v2_protein.faa.gz # Escherichia coli str. K-12 substr. MG1655 
// wget .../GCF_002803535.1_ASM280353v1_protein.faa.gz # Ochrobactrum pituitosum
// wget .../GCF_004208635.1_ASM420863v1_protein.faa.gz # Leucobacter triazinivorans
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/445/995/GCF_900445995.1_48290_B02/GCF_900445995.1_48290_B02_protein.faa.gz # Brevundimonas diminuta
// zcat dbs/custom/*.faa.gz | diamond makedb -p 8 -d custom --taxonmap prot.accession2taxid.gz --taxonnodes dbs/nodes.dmp

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()
reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/.ccs.fasta.gz/, file) }
assemblies = Channel.fromPath(params.assemblies, checkIfExists: true) 
                .map { file -> tuple(file.simpleName - ~/-hifiasm/, file) }

datasets = reads.join(assemblies)


process mask_assembly {
    tag "${strain}"

    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}.masked.fasta")

    script:
      """
      windowmasker -in $assembly \
                      -infmt fasta \
                      -mk_counts \
                      -sformat obinary \
                      -out tmp.counts 2> log \
        && windowmasker -in $assembly \
                        -infmt fasta \
                        -ustat tmp.counts \
                        -dust T \
                        -outfmt fasta \
                        -out ${strain}.masked.fasta
      """
}


process chunk_assembly {
    tag "${strain}"

    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}.chunks.fasta")

    script:
      """
      chunk_fasta.py --in ${assembly} \
        --chunk 100000 --overlap 0 --max-chunks 20 \
        --out ${strain}.chunks.fasta
      """
}

process diamond_search {
    tag "${strain}"

    input:
      tuple val(strain), path(assembly)
      path(dmnd_db)

    output:
      tuple val(strain), path("${strain}.chunk.diamond.tsv")

    script:
      """
      diamond blastx \
            --query ${assembly} \
            --db $params.dmnd_db \
            --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs $params.max_target_seqs \
            --max-hsps 1 \
            --evalue $params.evalue \
            --threads ${task.cpus} \
            > ${strain}.chunk.diamond.tsv
      """
}

process unchunk_hits {
    tag "${strain}"
    publishDir "$params.outdir/", mode: 'copy'

    input:
      tuple val(strain), path(diamondHits)

    output:
      tuple val(strain), path("${strain}.diamond.tsv")

    script:
      """
      unchunk_blast.py --in $diamondHits \
        --out ${strain}.diamond.tsv
      """
}

process map_reads {
    tag "${strain}"

    input:
      tuple val(strain), path(reads), path(assembly)

    output:
      tuple val(strain), path("${strain}.bam")

    script:
      """
      minimap2 -a -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 \
        --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 \
        --sam-hit-only -t ${task.cpus} ${assembly} \
        $reads | \
        samtools sort -@ ${task.cpus} -o ${strain}.bam
      """
}

process add_hits_and_coverage {
    tag "${strain}"
    publishDir "$params.outdir/btkDatasets", mode: 'copy'

    input:
      tuple val(strain), path(diamondHits), path(bam)

    output:
      tuple val(strain), path("${btkdir}")

    script:
      """
      $params.blobtoolsPath create \
            --fasta ${assembly} \
            --taxdump $params.taxdump \
            --taxid $params.taxid \
            $strain

      $params.blobtoolsPath add \
            --hits ${diamondHits} \
            --taxdump $params.taxdump \
            --taxrule $params.taxrule \
            --cov ${bam}=reads \
            $strain
      """
}


workflow {
    mask_assembly(assemblies) | chunk_assembly
    diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
    map_reads(datasets)
    add_hits_and_coverage(unchunk_hits.out.join(map_reads.out))
}

