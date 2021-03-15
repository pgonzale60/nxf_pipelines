nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "assemblyQC-${date}"
params.reads = "mini_DF5120.ccs.fasta.gz"
params.assemblies = "mini_DF5120.hifiasm.fasta.gz"
params.odb = 'nematoda_odb10'
params.busco_downloads = './busco_downloads'
params.telomere = 'TTAGGC'
params.busco2nigons = "gene2Nigon_busco20200927.tsv.gz"
params.min_occurr = 15
params.teloRepeatWindowSize = 1000
params.minimumGenesPerSequence = 15
params.minimumNigonFrac = 0.9
params.minimumFracAlignedTeloReads = 0.1
params.windowSizeQC = 5e5

reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(_filtered)?(\.telo)?(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }

fastFiles = Channel.fromPath(params.assemblies, checkIfExists: true)
assemblies = fastFiles.map { file -> tuple(file.Name - ~/(\.hifiasm)?(\.flye)?(\.wtdbg2)?(\.canu_plus_flye)?(\.canu)?(\.fa)?(\.fasta)?(\.gz)?$/, file.Name - ~/(\.fa)?(\.fasta)?(\.gz)?$/, file) }

busco2nigons = Channel.fromPath(params.busco2nigons, checkIfExists: true).collect()

busco_dbs = Channel.of(params.odb.split(','))
busco_db_dir = file(params.busco_downloads)
geno_busco = assemblies.combine(busco_dbs)



process busco {
    tag "${assembler}_${busco_db}"
    publishDir "$params.outdir/busco", mode: 'copy'

    input:
      tuple val(strain), val(assembler), path(genome), val(busco_db)
      path busco_db_dir

    output:
      path "*single_copy_busco_sequences.{faa,fna}"
      path "${assembler}_${busco_db}_short_summary.txt"
      tuple val(assembler), path( "${assembler}_${busco_db}_full_table.tsv"), emit: busco_full

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $genome > assembly.fasta
        else
            ln -s $genome assembly.fasta
      fi
      busco -c ${task.cpus} -l $busco_db -i assembly.fasta --out run_busco --mode geno
      awk 'BEGIN{FS="\\t";OFS=FS}(\$3 !~ /:/){print}' run_busco/run_*/full_table.tsv > ${assembler}_${busco_db}_full_table.tsv
      mv run_busco/short_summary* ${assembler}_${busco_db}_short_summary.txt
      #mv run_busco/run_*/full_table.tsv ${assembler}_${busco_db}_full_table.tsv
      for ext in .faa; do
        seqFile=${assembler}_${busco_db}_single_copy_busco_sequences\$ext
        for file in run_busco/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/*\$ext; do
          echo \">\$(basename \${file%\$ext})\" >> \$seqFile; tail -n +2 \$file >> \$seqFile;
        done
      done
      rm -rf run_busco/ assembly.fasta
      """
}

process get_telomeric_reads {
    tag "${strain}"
    label 'nemaQC'
    publishDir "$params.outdir/teloReads", mode: 'copy'

    input:
      tuple val(strain), path(reads)
    output:
      tuple val(strain), path("${strain}.telo.fasta.gz")

    script:
      """
      zcat $reads | filter_telomeric_reads.py --motif ${params.telomere} \
        --times ${params.min_occurr} --out ${strain}.telo.fasta.gz
      """
}

process map_telomeric_reads {
    tag "${assemblies[1]}"
    publishDir "$params.outdir/teloMaps", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(reads), val(assemblies)

    output:
      tuple val("${assemblies[1]}"), path("${assemblies[1]}.teloMapped.paf.gz")

    script:
      """
      minimap2 ${assemblies[2]} ${reads[1]} | \
        gzip -c > ${assemblies[1]}.teloMapped.paf.gz
      """
}


process count_telomeric_repeat {
    tag "${assembler}"
    publishDir "$params.outdir/teloRepeatCounts", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(strain), val(assembler), path(assembly)
    
    output:
      tuple val(assembler), path( "${assembler}_teloRepeatCounts.tsv.gz")

    script:
      """
      zcat $assembly | \
        awk '/^>/{if (l!="") print l; printf "%s\\t", \$1; l=0; next}{l+=length(\$0)}END{print l}' | \
        sed "s/>//" > ${assembler}.seqlen.tsv
      seqkit locate --bed -M -G -p ${params.telomere} $assembly | \
        cut -f 1,2,3 | \
        sort -k1,1 -k2,2n > ${assembler}.teloRepeats.tsv
      bedtools makewindows -g ${assembler}.seqlen.tsv -w ${params.teloRepeatWindowSize} | \
        bedtools intersect -a stdin -b ${assembler}.teloRepeats.tsv -wa -wb | \
        bedtools groupby -i stdin -g 1,2,3 -c 1 -o count | \
        awk -F '\\t' 'BEGIN{OFS=FS}{print \$0, "telomeres"}' | \
        gzip -c > ${assembler}_teloRepeatCounts.tsv.gz
      rm ${assembler}.seqlen.tsv ${assembler}.teloRepeats.tsv
      """
}


process nematode_chromosome_QC {
    tag "${assembler}"
    publishDir "$params.outdir/nemaChromQC", mode: 'copy'
    label 'nemaQC'

    input:
      tuple val(assembler), path(buscoTable), path(teloMappedReads), path(teloRepeats)
      path(busco2nigons)
    
    output:
      path "${assembler}.{pdf,buscoString.txt,teloMappedBlocks.tsv}"
      path "${assembler}.chromQC.tsv" , emit: busco_full

    script:
      """
      nemaChromQC.R --assemblyName $assembler \
        --nigon $busco2nigons --busco $buscoTable \
        --teloMappedPaf $teloMappedReads \
        --teloRepeats $teloRepeats \
        --minimumGenesPerSequence $params.minimumGenesPerSequence \
        --minimumNigonFrac $params.minimumNigonFrac \
        --minimumFracAlignedTeloReads $params.minimumFracAlignedTeloReads \
        --windowSize $params.windowSizeQC
      """
}

process get_contiguity_stats {
    tag "all"
    publishDir "$params.outdir/", mode: 'copy'
    label 'nemaQC'

    input:
      path(assemblies)
    
    output:
      path "assemblies_contiguity_stats.tsv"

    script:
      """
      seqkit stats -a -T $assemblies > assemblies_contiguity_stats.tsv
      """
}


workflow {
    busco(geno_busco, busco_db_dir)
    count_telomeric_repeat(assemblies)
    get_telomeric_reads(reads)
    map_telomeric_reads(get_telomeric_reads.out.cross(assemblies))
    get_contiguity_stats(fastFiles.collect())
    nematode_chromosome_QC(busco.out.busco_full.join(map_telomeric_reads.out.join(count_telomeric_repeat.out)),
     busco2nigons)
}
