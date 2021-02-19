nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "teloAssem-${date}"
params.telomere = 'TTAGGC'
params.min_occurr = 15
params.reads = "mini_DF5120.ccs.fasta.gz"


reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }



process get_telomeric_reads {
    tag "${strain}"
    label 'btk'
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


process hifiasm {
    tag "${strain}"
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

process direct_by_telomere {
    tag "${strain}"
    publishDir "$params.outdir/assemblies", mode: 'copy'

    input:
      tuple val(strain), path(assembly)

    output:
      tuple val(strain), path("${strain}.hifiasm.telopointed.fasta.gz")

    script:
      """
      seqkit locate --bed -M -G -p ${params.telomere} $assembly > tmp
      bioawk -t '{tlen[\$1]+=\$2; nr[\$1]+=1}END { for (i in tlen) print i, tlen[i]/nr[i]}' tmp | sort -k1,1 > avgTeloOccur.tsv
      zcat $assembly | awk '/^>/{if (l!="") print l; printf "%s\\t", \$1; l=0; next}{l+=length(\$0)}END{print l}' | sed "s/>//" | sort -k1,1 > seqlen.tsv
      join seqlen.tsv avgTeloOccur.tsv > both.tsv
      awk '\$3>\$2/2{print \$1}' both.tsv > idsToRev.txt
      awk '\$3<=\$2/2{print \$1}' both.tsv > straightIds.txt
      seqtk subseq $assembly idsToRev.txt | seqkit seq -r | gzip -c > reversed.fasta.gz
      seqtk subseq -l 60 $assembly straightIds.txt | gzip -c > straight.fasta.gz
      cat reversed.fasta.gz straight.fasta.gz > ${strain}.hifiasm.telopointed.fasta.gz
      """
}

process map_reads {
    tag "${strain}"
    publishDir "$params.outdir/bams", mode: 'copy'
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

process coverage_dist {
    tag "${strain}"
    publishDir "$params.outdir/cov_depth", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(bam)

    output:
      tuple val(strain), path("${strain}.depth.tsv.gz")

    script:
      """
      samtools depth $bam | gzip -c > ${strain}.depth.tsv.gz
      """
}

workflow {
    get_telomeric_reads(reads) //| hifiasm | direct_by_telomere
    // map_reads(get_telomeric_reads.out.join(direct_by_telomere.out)) | coverage_dist
}

/*
cut -f 1 telomeres_flye_CEW1.fasta.fai > ids
printf "" > allout

while read p; do
samtools faidx telomeres_flye_CEW1.fasta $p > query
grep -v $p ids | xargs samtools faidx telomeres_flye_CEW1.fasta > target
minimap2 target query >> allout
echo $p
done < ids

bioawk -t '{print $1, $3, $4, $6, $8, $9}' allout > data/selfmap.txt

bioawk -t '{print "chr", "-", $1, $1, 0, $2, "lblue"}' telomeres_flye_CEW1.fasta.fai > data/otitelo_karyotype.txt


## gscope parsing #
printf "%s\tkcov\t" sample > miniBtk-20210123/gscope.tsv
head -n 1 miniBtk-20210123/genomescope/PS2068_k31_gscope/summary.tsv | tr '\n' '\t' >> miniBtk-20210123/gscope.tsv # tr is to remove end of line from and keep writing to the header on the next line
printf "plot\tlog plot\n" ${sampName} ${sampName} >> miniBtk-20210123/gscope.tsv
for fdir in miniBtk-20210123/genomescope/*; do
sampName=$(basename $fdir)
printf "%s\t" ${sampName%_k31_gscope} >> miniBtk-20210123/gscope.tsv
grep "^kmercov" ${fdir}/model.txt | tr -s " " | cut -f 2 -d ' ' |  xargs printf "%.1f\t" >> miniBtk-20210123/gscope.tsv
tail -n 1 ${fdir}/summary.tsv | tr '\n' '\t' >> miniBtk-20210123/gscope.tsv
printf "assets/images/kmers/%s/plot.png\tassets/images/kmers/%s/plot.log.png\n" ${sampName} ${sampName} >> miniBtk-20210123/gscope.tsv
done

## Btk weighted average coverage #


printf "sample\tweighted_average_coverage\n" > bothBtks_cov.tsv
for dataset in ~/testNext/bothBtks/*filtered; do
sampName=$(basename  $dataset)
~/sw/blobtoolkit/blobtools2/blobtools filter --table exe.tsv  $dataset
printf "%s\t" $sampName >> bothBtks_cov.tsv
bioawk -t 'BEGIN{tsum=0; tlen=0}NR>1{tsum+=$5*$4;tlen+=$4}END{printf "%.1f\n", tsum/tlen}' exe.tsv >> bothBtks_cov.tsv
done

*/
