nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "miniBtk-${date}"
params.reads = "/home/ubuntu/oscheius/0-inputs/test.ccsf.fasta.gz"
params.assemblies = "/home/ubuntu/oscheius/0-inputs/test.fasta"
params.odb = 'nematoda_odb10'
params.busco_downloads = "/lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/dbs/busco_downloads/"
params.dmnd_db = "custom.dmnd"
params.sampRate = 0.2
params.max_target_seqs = 1000
params.evalue = 0.000001
params.taxid = 2613844
params.taxdump = "taxdump/"
params.taxrule = "bestsumorder"




// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/067/145/GCA_904067145.1_BOKI2/GCA_904067145.1_BOKI2_protein.faa.gz # Bursaphelenchus okinawaensis
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz # Caenorhabditis elegans
// wget .../GCF_000005845.2_ASM584v2_protein.faa.gz # Escherichia coli str. K-12 substr. MG1655 
// wget .../GCF_002803535.1_ASM280353v1_protein.faa.gz # Ochrobactrum pituitosum
// wget .../GCF_004208635.1_ASM420863v1_protein.faa.gz # Leucobacter triazinivorans
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/445/995/GCF_900445995.1_48290_B02/GCF_900445995.1_48290_B02_protein.faa.gz # Brevundimonas diminuta
// zcat dbs/custom/*.faa.gz | diamond makedb -p 8 -d custom --taxonmap prot.accession2taxid.gz --taxonnodes dbs/nodes.dmp

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()
taxdump_db = Channel.fromPath(params.taxdump, checkIfExists: true).collect()
reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_hifi_reads.fasta.gz)?(\.merged)?(\.subsamp)?(\.ccs)?(\.fa)?(\.fasta)?(\.fastq)?(\.gz)?$/, file) }
assemblies = Channel.fromPath(params.assemblies, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.fasta)?(\.gz)?$/, file.Name - ~/(_hifiasm.bp.p_ctg.fa)?(\.noTelos)?(\.flyemeta)?(\.hifiasm)?(\.tol)?(\.g30k)?(\.purged)?(\.l500k)?(\.salsa)?(\.juiced)?(\.noCont)?(\.noMito)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }
busco_dbs = Channel.of(params.odb.split(','))
busco_db_dir = file(params.busco_downloads)


process busco {
    tag "${assembler}_${busco_db}"
    publishDir "$params.outdir/busco", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), val(assembler), path(genome), val(busco_db)
      path busco_db_dir

    output:
      path "*single_copy_busco_sequences.{faa,fna}"
      path "${assembler}_${busco_db}_short_summary.txt"
      tuple val(assembler), path( "${assembler}_${busco_db}_full_table.tsv"), emit: busco_full

    script:
      """
      export http_proxy=http://wwwcache.sanger.ac.uk:3128
      export https_proxy=http://wwwcache.sanger.ac.uk:3128
      busco -c ${task.cpus} -l $busco_db -i $genome --out run_busco --mode geno
      awk 'BEGIN{FS="\\t";OFS=FS}(\$3 !~ /:/){print}' run_busco/run_*/full_table.tsv > ${assembler}_${busco_db}_full_table.tsv
      mv run_busco/short_summary* ${assembler}_${busco_db}_short_summary.txt
      #mv run_busco/run_*/full_table.tsv ${assembler}_${busco_db}_full_table.tsv
      for ext in .faa; do
        seqFile=${assembler}_${busco_db}_single_copy_busco_sequences\$ext
        for file in run_busco/run_${busco_db}/busco_sequences/single_copy_busco_sequences/*\$ext; do
          echo \">\$(basename \${file%\$ext})\" >> \$seqFile; tail -n +2 \$file >> \$seqFile;
        done
      done
      rm -rf run_busco/
      """
}

process mask_assembly {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), val(assemName), path(assembly)

    output:
      tuple val(strain), val(assemName), path("${strain}.masked.fasta")

    script:
      """
      cp $assembly assembly.fasta
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
                        -out ${strain}.masked.fasta
      rm assembly.fasta
      """
}


process chunk_assembly {
    tag "${assemName}"
    label 'btk'

    input:
      tuple val(assemName), val(strain), path(assembly)

    output:
      tuple val(assemName),path("${strain}.chunks.fasta")

    script:
      """
      chunk_fasta.py --in ${assembly} \
        --chunk 100000 --overlap 0 --max-chunks 20 \
        --out ${strain}.chunks.fasta
      """
}

process diamond_search {
    tag "${strain}"
    label 'big_parallelizable'
    label 'btk'

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
    label 'btk'

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
    tag "${assemblies[1]}"
    publishDir "$params.outdir/", mode: 'copy'

    input:
      tuple val(reads), val(assemblies)

    output:
      tuple val("${assemblies[1]}"), path("${assemblies[1]}.bam")

    script:
      """
      /software/team301/minimap2-2.24/minimap2 -ax map-hifi \
        --sam-hit-only -t ${task.cpus} ${assemblies[2]} \
        ${reads[1]} | \
        samtools sort -@ ${task.cpus} -o ${assemblies[1]}.bam
      """
}

process add_hits_and_coverage {
    tag "${strain}"
    publishDir "$params.outdir/btkDatasets", mode: 'copy'
    label 'btk'

    input:
      tuple val(assemName), val(strain), path(assembly), path(diamondHits), path(bam), path(busco_full_tsv)
      path(taxdump_db)

    output:
      tuple val(assemName), path("${assemName}"), emit: blobDir

    script:
      """
      blobtools create \
            --fasta ${assembly} \
            --taxdump $taxdump_db \
            --taxid $params.taxid \
            $assemName

      blobtools add \
            --hits ${diamondHits} \
            --busco ${busco_full_tsv} \
            --taxdump $taxdump_db \
            --taxrule $params.taxrule \
            --cov ${bam} \
            $assemName
      """
}


workflow {
    mask_assembly(assemblies) | chunk_assembly
    busco(assemblies.combine(busco_dbs), busco_db_dir)
    diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
    map_reads(reads.cross(assemblies.map{it -> tuple(it[1], it[0], it[2]) }))
    add_hits_and_coverage(assemblies.join(unchunk_hits.out.join(map_reads.out.join(busco.out.busco_full))), taxdump_db) 
}
