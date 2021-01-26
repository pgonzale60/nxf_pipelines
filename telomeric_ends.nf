nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "teloAssem-${date}"
params.telomere = 'TTAGGCTTAGGCTTAGGCTTAGGCTTAGGC'
params.min_occurr = 1
params.reads = "mini_DF5120.ccs.fasta.gz"


reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }



process get_telomeric_reads {
    tag "${strain}"

    input:
      tuple val(strain), path(reads)
    output:
      tuple val(strain), path("${strain}.telo.fasta.gz")

    script:
      """
      seqkit locate --bed -M -G -p ${params.telomere} $reads | \
      cut -f 1 | uniq -c | tr -s " " | awk '\$1>${params.min_occurr}{print \$2}' > teloReads.ids
      seqtk subseq $reads teloReads.ids | gzip -c > ${strain}.telo.fasta.gz
      """
}


process hifiasm {
    tag "${strain}"
    publishDir "$params.outdir/assemblies", mode: 'copy'
    label 'big_parallelizable'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.hifiasm.fasta.gz")

    script:
      """
      /software/team301/hifiasm/hifiasm $reads -o $strain -t ${task.cpus}
      awk '/^S/{print ">"\$2"\\n"\$3}' ${strain}.p_ctg.gfa | fold | gzip -c > ${strain}.hifiasm.fasta.gz
      """
}

workflow {
    get_telomeric_reads(reads) | hifiasm
}