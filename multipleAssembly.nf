nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "multiassem-${date}"
params.reads = "mini_DF5120.ccs.fasta.gz"


reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }

/*
assemblies = Channel
                .fromPath(params.assemblies, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.fa)?(\.fasta)?(\.gz)?$/, file) }

dmnd_db = Channel.fromPath(params.dmnd_db, checkIfExists: true).collect()
*/


process hifiasm {
    tag "${strain}"
    publishDir "$params.outdir/hifiasm", mode: 'copy'

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

process flye {
    tag "${strain}"
    publishDir "$params.outdir", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.flye.fasta.gz")

    script:
      """
      /software/team301/Flye-2.8.2/Flye/bin/flye --threads ${task.cpus} \
       --pacbio-hifi $reads --meta -o flyemeta
       cat flyemeta/assembly.fasta | bgzip -c > ${strain}.flye.fasta.gz
      """
}

process wtdbg2 {
    tag "${strain}"
    publishDir "$params.outdir", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.wtdbg2.fasta.gz")

    script:
      """
      /software/team301/wtdbg2/wtdbg2.pl -t ${task.cpus} -x ccs -g 100m -o ${strain} $reads
       bgzip ${strain}.cns.fa
       mv ${strain}.cns.fa.gz ${strain}.wtdbg2.fasta.gz
      """
}

process canu {
    tag "${strain}"
    publishDir "$params.outdir", mode: 'symlink'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}")

    script:
      """
      /software/tola/bin/canu-2.1.1/bin/canu -d ${strain} \
        -p ${strain} gridEngineResourceOption='-M MEMORY -R "select[mem>MEMORY] rusage[mem=MEMORY]" -n THREADS -R "span[hosts=1]"' \
        genomeSize=60M -pacbio-hifi $reads
      """
}





workflow {
    flye(reads)
    wtdbg2(reads)
    canu(reads)
    /*mask_assembly(assemblies) | chunk_assembly
    diamond_search(chunk_assembly.out, dmnd_db) | unchunk_hits
    */
}
