nextflow.preview.dsl=2

// Selecting nematoda (6231) genomes from NCBI with dataset v6.3.0
// datasets assembly-descriptors tax-id 6231 > tm.json
// assemblies with .assembly_category != null are either representative or reference genomes
// jq -r '.datasets[] | select(.assembly_category != null) | [.assembly_accession, .org.sci_name] | @tsv' tm.json | sort -u -k2,2 -t $'\t' | sed 's/sp./sp/; s/ /_/g' > selected_assems.tsv

params.accessions = 'accessions.txt'
params.odb = 'nematoda_odb10'
params.busco_downloads = './busco_downloads'
params.outdir = './results'



busco_db = params.odb

busco_db_dir = Channel.fromPath(params.busco_downloads, type: 'dir', checkIfExists: true).collect()

accessions = Channel
   .fromPath("${params.accessions}", checkIfExists: true)
   .splitCsv( by: 1, sep: '\t' )



process downAssem {
   tag "${sci_name}"
   publishDir "$params.outdir/assemblies", mode: 'copy'

    input:
         tuple val(assemLink), val(sci_name)

    output:
        tuple val(sci_name), path("${sci_name}.fasta.gz"), emit: input_fastas

    script:
    """
        curl ${assemLink} > ${sci_name}.fasta.gz
    """

}

process busco {
    tag "${sci_name}"
    publishDir "$params.outdir/${sci_name}", mode: 'copy'

    input:
      tuple val(sci_name), path(genome)
      val busco_db
      path busco_db_dir

    output:
      path "${sci_name}_${busco_db}_full_table.tsv", emit: full_busco_table
      path "${sci_name}_${busco_db}_short_summary.txt", emit: short_busco_report
      path "${sci_name}_${busco_db}_single_copy_busco_sequences*", emit: busco_busco_sequences

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $genome > assembly.fasta
        else
            ln -s $genome assembly.fasta
      fi
      busco -c ${task.cpus} -l $busco_db -i assembly.fasta --out run_busco --mode geno
      mv run_busco/short_summary* ${sci_name}_${busco_db}_short_summary.txt
      awk 'BEGIN{FS="\\t";OFS=FS}(\$3 !~ /:/){print}' run_busco/run_*/full_table.tsv > ${sci_name}_${busco_db}_full_table.tsv
      #mv run_busco/run_*/full_table.tsv ${sci_name}_${busco_db}_full_table.tsv

      for ext in .faa; do
        seqFile=${sci_name}_${busco_db}_single_copy_busco_sequences\$ext
        for file in run_busco/run_${busco_db}/busco_sequences/single_copy_busco_sequences/*\$ext; do \
          echo \">\$(basename \${file%\$ext})\" >> \$seqFile; tail -n +2 \$file >> \$seqFile;
        done
      done
      rm -rf run_busco/run_${busco_db}/ run_busco/blast_db run_busco/run_*/blast_output run_busco/run_*/hmmer_output assembly.fasta
      """
}


workflow {
    downAssem(accessions)
    busco(downAssem.out, busco_db, busco_db_dir)
}
