nextflow.enable.dsl=2



params.odb = 'lepidoptera_odb10,nematoda_odb10'
params.busco_downloads = './busco_downloads'
params.genomes = './data/*fasta.gz'
params.outdir = './results'



busco_dbs = Channel.of(params.odb.split(','))

busco_db_dir = file(params.busco_downloads)

genomes = Channel
                .fromPath(params.genomes, checkIfExists: true)
                .map { file -> tuple(file.simpleName, file) }




geno_busco = genomes.combine(busco_dbs)



process busco {
    tag "${assembler}_${busco_db}"
    publishDir "$params.outdir/busco", mode: 'copy'

    input:
      tuple val(assembler), path(genome), val(busco_db)
      path busco_db_dir

    output:
      path "*single_copy_busco_sequences.faa"
      path "${assembler}_${busco_db}_short_summary.txt"
      tuple val(assembler), path( "${assembler}_${busco_db}_full_table.tsv"), emit: busco_full

    script:
      """
      export http_proxy=http://wwwcache.sanger.ac.uk:3128
      export https_proxy=http://wwwcache.sanger.ac.uk:3128
      if [ -f *.gz ]; then
            gunzip -c $genome > assembly.fasta
        else
            ln $genome assembly.fasta
      fi
      busco -c ${task.cpus} -l $busco_db -i assembly.fasta --out run_busco --mode geno --offline
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

workflow {
    busco(geno_busco, busco_db_dir)
}