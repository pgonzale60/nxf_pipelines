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
params.kmer = '31'
params.odb = 'nematoda_odb10'
params.busco_downloads = '/lustre/scratch116/tol/teams/team301/dbs/busco_2020_08/busco_downloads/'



// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/904/067/145/GCA_904067145.1_BOKI2/GCA_904067145.1_BOKI2_protein.faa.gz # Bursaphelenchus okinawaensis
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz # Caenorhabditis elegans
// wget .../GCF_000005845.2_ASM584v2_protein.faa.gz # Escherichia coli str. K-12 substr. MG1655 
// wget .../GCF_002803535.1_ASM280353v1_protein.faa.gz # Ochrobactrum pituitosum
// wget .../GCF_004208635.1_ASM420863v1_protein.faa.gz # Leucobacter triazinivorans
// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/445/995/GCF_900445995.1_48290_B02/GCF_900445995.1_48290_B02_protein.faa.gz # Brevundimonas diminuta
// zcat dbs/custom/*.faa.gz | diamond makedb -p 8 -d custom --taxonmap prot.accession2taxid.gz --taxonnodes dbs/nodes.dmp

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()
busco_db_dir = Channel.fromPath(params.busco_downloads, checkIfExists: true, type: 'dir').collect()
busco_dbs = Channel.of(params.odb.split(','))

reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }



process kmer_hist {
    tag "${strain}"
    publishDir "$params.outdir/kat", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads)
    output:
      tuple val(strain), path("${strain}.hist"), emit: kmer_counts
      path("${strain}.hist*{json,png}")


    script:
      """
      kat hist -o ${strain}.hist -t ${task.cpus} -m ${params.kmer} $reads
      """
}


process genomescope {
    tag "${strain}"
    publishDir "$params.outdir/genomescope", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(histo)
    output:
      path("${strain}_k${params.kmer}_gscope")
      
    script:
      """
      genomescope.R $histo $params.kmer 150 ${strain}_k${params.kmer}_gscope
      """
}


process hifiasm {
    tag "${strain}"
    publishDir "$params.outdir/assemblies", mode: 'copy'
    label 'big_parallelizable'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.hifiasm.fasta")

    script:
      """
      /software/team301/hifiasm/hifiasm $reads -o $strain -t ${task.cpus}
      awk '/^S/{print ">"\$2"\\n"\$3}' ${strain}.p_ctg.gfa | fold > ${strain}.hifiasm.fasta
      """
}

process busco {
    tag "${strain}_${busco_db}"
    publishDir "$params.outdir/busco/${strain}_${busco_db}", mode: 'copy'
    label 'big_parallelizable'

    input:
      tuple val(strain), path(genome), val(busco_db)
      path busco_db_dir

    output:
      path("${strain}_${busco_db}_single_copy_busco_sequences*")
      tuple val(strain), path("${strain}_${busco_db}_full_table.tsv"), emit: busco_table
      path("${strain}_${busco_db}_short_summary.txt")

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $genome > assembly.fasta
        else
            ln $genome assembly.fasta
      fi
      export AUGUSTUS_CONFIG_PATH=augustus_conf
      cp -r /augustus/config/ \$AUGUSTUS_CONFIG_PATH
      busco -c ${task.cpus} -l $busco_db -i assembly.fasta --out run_busco --mode geno
      mv run_busco/short_summary* ${strain}_${busco_db}_short_summary.txt
      mv run_busco/run_*/full_table.tsv ${strain}_${busco_db}_full_table.tsv
      for ext in .faa .fna; do
        seqFile=${strain}_${busco_db}_single_copy_busco_sequences\$ext
        for file in run_busco/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/*\$ext; do
          echo \">\$(basename \${file%\$ext})\" >> \$seqFile; tail -n +2 \$file >> \$seqFile;
        done
      done
      rm -rf \$AUGUSTUS_CONFIG_PATH run_busco/ assembly.fasta
      """
}

process mask_assembly {
    tag "${strain}"
    label 'btk'

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
    label 'btk'

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
    tag "${strain}"
    label 'big_parallelizable'
    label 'btk'

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

process kat_plot {
    tag "${strain}"
    publishDir "$params.outdir/kat", mode: 'copy'
    label 'big_parallelizable'
    label 'btk'

    input:
      tuple val(strain), path(reads), path(assembly)

    output:
      tuple val(strain), path("${strain}-main.mx.*")

    script:
      """
      kat comp -t ${task.cpus} -o $strain $reads $assembly
      """
}

process create_blobDir {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}")

    script:
      """
echo \$PATH
      $params.blobtoolsPath create \
            --fasta ${assembly} \
            --taxdump $params.taxdump \
            --taxid $params.taxid \
            $strain
      """
}

process add_hits_coverage_and_busco {
    tag "${strain}"
    publishDir "$params.outdir/btkDatasets", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(btkdir), path(diamondHits), path(bam), path(busco)

    output:
      tuple val(strain), path("${btkdir}")

    script:
      """
      $params.blobtoolsPath add \
            --hits $diamondHits \
            --taxdump $params.taxdump \
            --taxrule $params.taxrule \
            --cov ${bam}=reads \
            --busco $busco \
            $btkdir

      $params.blobtoolsPath filter \
            --summary $btkdir/summary.json \
            $btkdir
      """
}

process btk_static_images {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(btkdir)

    output:
      tuple val(strain), path("${btkdir}")

    script:
      """
      $params.blobtoolsPath view \
            --view blob \
            --param plotShape=circle \
            --param bestsumorder_phylum--Order=no-hit%2CNematoda%2CProteobacteria%2CActinobacteria \
            --format png --format svg \
            $btkdir
      $params.blobtoolsPath view \
            --view cumulative \
            --format png --format svg \
            $btkdir
      mv *svg *png $btkdir
      """
}

process filter_fasta {
    tag "${strain}"
    publishDir "$params.outdir/filteredData", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(btkdir), path(bam), path(assembly), path(reads)

    output:
      tuple val(strainName), path("$filtered_assemFile"), emit: filtered_assem
      tuple val(strainName), path("$filtered_readsFile"), emit: filtered_reads

    script:
    strainName = strain + "_filtered"
    btk_fltrd_assemFile = assembly.baseName - ~/(\.fasta)?(\.fa)?$/ + ".filtered.fasta"
    btk_fltrd_readsFile = reads.baseName - ~/(\.gz)?$/ + ".filtered.gz"
    filtered_assemFile = strainName + ".hifiasm.fasta"
    filtered_readsFile = strainName + ".ccs.fasta.gz"
      """
      $params.blobtoolsPath filter \
        --param bestsumorder_superkingdom--Keys=Bacteria \
        --param gc--Max=0.49 \
        --fasta $assembly \
        --fastq $reads \
        --cov $bam \
        $btkdir
      mv $btk_fltrd_readsFile $filtered_readsFile
      mv $btk_fltrd_assemFile $filtered_assemFile
      """
}


workflow raw_asses {
    take: reads
    main:
        kmer_hist(reads)
        genomescope(kmer_hist.out.kmer_counts)
        hifiasm(reads) | mask_assembly | chunk_assembly
        busco(hifiasm.out.combine(busco_dbs), busco_db_dir)
        diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
        map_reads(reads.join(hifiasm.out))
        kat_plot(reads.join(hifiasm.out))
        create_blobDir(hifiasm.out)
        add_hits_coverage_and_busco(create_blobDir.out.join(unchunk_hits.out.join(map_reads.out.join(busco.out.busco_table))))
        filter_fasta(add_hits_coverage_and_busco.out.join(map_reads.out.join(hifiasm.out.join(reads))))
    emit:
        filter_fasta.out.filtered_reads
}

workflow fltd_asses {
    take: reads
    main:
        kmer_hist(reads)
        genomescope(kmer_hist.out.kmer_counts)
        hifiasm(reads) | mask_assembly | chunk_assembly
        busco(hifiasm.out.combine(busco_dbs), busco_db_dir)
        diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
        map_reads(reads.join(hifiasm.out))
        kat_plot(reads.join(hifiasm.out))
        create_blobDir(hifiasm.out)
        add_hits_coverage_and_busco(create_blobDir.out.join(unchunk_hits.out.join(map_reads.out.join(busco.out.busco_table))))
    emit:
        add_hits_coverage_and_busco.out
}

workflow {
    raw_asses(reads)
    fltd_asses(raw_asses.out)
}