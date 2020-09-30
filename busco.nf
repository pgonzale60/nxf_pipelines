nextflow.enable.dsl=2



params.odb = 'lepidoptera_odb10,endopterygota_odb10'
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
    tag "${species}_${busco_db}"
    publishDir params.outdir

    container = "$HOME/busco_v4.1.2_cv1.sif"
    cpus 8

    input:
      tuple val(species), path(genome), val(busco_db)
      path busco_db_dir
    output:
      path "${species}_full_table.tsv", emit: full_busco_table
      path "short_summary_${species}.txt", emit: short_busco_report

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
      mv run_busco/short_summary* short_summary_${species}.txt
      mv run_busco/run_*/full_table.tsv ${species}_full_table.tsv
      rm -rf \$AUGUSTUS_CONFIG_PATH run_busco/* assembly.fasta
      """
}

workflow {
    busco(geno_busco, busco_db_dir)
}