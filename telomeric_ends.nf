nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "teloAssem-${date}"
params.telomere = 'TTAGGC'
params.min_occurr = 15
params.reads = "mini_DF5120.ccs.fasta.gz"
params.memeMotif = "DF5120_diminution_400.meme.txt"
params.teloRepeatWindowSize = 1000
params.dimiWidth = 200


reads = Channel.fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.Name - ~/(\.ccs)?(\.fa)?(\.fasta)?(\.gz)?$/, file) }

memeMotif = Channel.fromPath(params.memeMotif, checkIfExists: true).collect()

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
    publishDir "$params.outdir/dimiAssem", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads)

    output:
      tuple val(strain), path("${strain}.hifiasm.fasta.gz")

    script:
      """
      /software/team301/hifiasm/hifiasm $reads -o $strain -t ${task.cpus}
      awk '/^S/{print ">"\$2"\\n"\$3}' ${strain}.p_ctg.gfa | fold | bgzip -c > ${strain}.hifiasm.fasta.gz
      """
}

process flye {
    tag "${strain}"
    publishDir "$params.outdir/dimiAssem", mode: 'copy'
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

process get_read_coordinates {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(bam)

    output:
      tuple val(strain), path("${strain}.teloMapped.coords.tsv")

    script:
      """
      samtools view ${bam} | teloCoord > ${strain}.teloMapped.coords.tsv
      """
}

process group_read_coordinates {
    tag "${strain}"
    label 'r'
    publishDir "$params.outdir/telomere_read_coords", mode: 'copy'

    input:
      tuple val(strain), path(coords)

    output:
      tuple val(strain), path("${strain}.teloPositions.tsv")

    script:
      """
      grCoords.R ${coords} ${strain}.teloPositions.tsv
      """
}

process get_diminuted_regions {
    tag "${strain}"
    label 'btk'

    input:
      tuple val(strain), path(teloPos), path(fasta)

    output:
      tuple val(strain), path("${strain}.diminuted.fasta")

    script:
      """
      awk 'BEGIN{FS="\\t"}(NR>1){startCoord=\$2-${params.dimiWidth}; if(0 > startCoord){startCoord=0}; print \$1":"startCoord"-"\$2+${params.dimiWidth}}' $teloPos | \
        xargs samtools faidx $fasta > ${strain}.diminuted.fasta
      """
}

process map_telomeric_reads {
    tag "$strain"
    publishDir "$params.outdir/teloMaps", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(reads), path(assembly)

    output:
      tuple val("$strain"), path("${strain}.teloMapped.paf.gz")

    script:
      """
      minimap2 $assembly $reads | \
        gzip -c > ${strain}.teloMapped.paf.gz
      """
}

process bam2fasta {
    tag "$strain"
    label 'btk'

    input:
      tuple val(strain), path(bam)

    output:
      tuple val("$strain"), path("${strain}.mappedReads.fasta.gz")

    script:
      """
      samtools fasta $bam | bgzip -c > ${strain}.mappedReads.fasta.gz
      """  
}

process bam2coords {
    tag "$strain"

    input:
      tuple val(strain), path(bam)

    output:
      tuple val("$strain"), path("${strain}.mappedCoords.tsv.gz")

    script:
      """
      samtools view $bam | teloCoord | bgzip -c > ${strain}.mappedCoords.tsv.gz
      """  
}

process remove_end_reads {
    tag "$strain"
    label 'btk'

    input:
      tuple val(strain), path(reads), path(end_reads)

    output:
      tuple val("$strain"), path("${strain}.dimiReads.fasta.gz")

    script:
      """
      seqkit seq $end_reads -n > end_reads.ids.txt
      seqkit seq $reads -n | sort > some_reads.ids.txt
      cat end_reads.ids.txt some_reads.ids.txt | sort | \
        uniq --unique > unique_reads.txt
      # We need to join becuase there can be end_reads that did not map
      # to the assembly and hence will appear only once in the above list
      join some_reads.ids.txt unique_reads.txt > non_end_reads.txt
      cat non_end_reads.txt | xargs samtools faidx $reads | \
        bgzip -c > ${strain}.dimiReads.fasta.gz
      """  
}

process fimo {
    tag "$strain"
    publishDir "$params.outdir/fimo", mode: 'copy'

    input:
      tuple val(strain), path(assembly)
      path(motif)

    output:
      tuple val("$strain"), path("${strain}.fimo.tsv")

    script:
      """
      if [ -f *.gz ]; then
            gunzip -c $assembly > assembly.fasta
        else
            ln -s $assembly assembly.fasta
      fi
      fimo --max-strand $motif assembly.fasta
      mv fimo_out/fimo.tsv ${strain}.fimo.tsv
      """  
}

process meme {
    tag "$strain"
    publishDir "$params.outdir/meme", mode: 'copy'

    input:
      tuple val(strain), path(diminuted_regions)

    output:
      path "${strain}.meme.{txt,html}"

    script:
      """
      meme -dna $diminuted_regions
      mv meme_out/meme.txt ${strain}.meme.txt
      mv meme_out/meme.html ${strain}.meme.html
      """  
}

process count_telomeric_repeat {
    tag "${strain}"
    publishDir "$params.outdir/teloRepeatCounts", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(assembly)
    
    output:
      tuple val(strain), path( "${strain}_teloRepeatCounts.tsv.gz")

    script:
      """
      zcat $assembly | \
        awk '/^>/{if (l!="") print l; printf "%s\\t", \$1; l=0; next}{l+=length(\$0)}END{print l}' | \
        sed "s/>//" > ${strain}.seqlen.tsv
      seqkit locate --bed -M -G -p ${params.telomere} $assembly | \
        cut -f 1,2,3 | \
        sort -k1,1 -k2,2n > ${strain}.teloRepeats.tsv
      bedtools makewindows -g ${strain}.seqlen.tsv -w ${params.teloRepeatWindowSize} | \
        bedtools intersect -a stdin -b ${strain}.teloRepeats.tsv -wa -wb | \
        bedtools groupby -i stdin -g 1,2,3 -c 1 -o count | \
        awk -F '\\t' 'BEGIN{OFS=FS}{print \$0, "telomeres"}' | \
        gzip -c > ${strain}_teloRepeatCounts.tsv.gz
      rm ${strain}.seqlen.tsv ${strain}.teloRepeats.tsv
      """
}

process crossCheckMotif {
    tag "${strain}"
    publishDir "$params.outdir/crossCheckMotif", mode: 'copy'
    label 'btk'

    input:
      tuple val(strain), path(mappedReads), path(motifCoords)
    
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

workflow get_diminuted_reads {
  take:
    reads
  main:
    get_telomeric_reads(reads) | hifiasm //| direct_by_telomere
    map_reads(reads.join(hifiasm.out)) | bam2fasta
    remove_end_reads(bam2fasta.out.join(get_telomeric_reads.out))
  emit: 
    diminuted_reads = remove_end_reads.out
    end_reads = get_telomeric_reads.out
    end_related_reads = bam2fasta.out
}

workflow locate_diminution {
  take:
    diminuted_reads
    end_reads
    memeMotif
    end_related_reads
  main:
    flye(end_related_reads) | count_telomeric_repeat
    fimo(flye.out, memeMotif)
    map_reads(end_reads.join(flye.out)) | get_read_coordinates | group_read_coordinates
    get_diminuted_regions(group_read_coordinates.out.join(flye.out)) | meme
  emit: 
    meme.out
}

workflow {
    get_diminuted_reads(reads)
    locate_diminution(get_diminuted_reads.out.diminuted_reads, \
      get_diminuted_reads.out.end_reads, \
      memeMotif, get_diminuted_reads.out.end_related_reads)
}

// 

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


nextflow -C /lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/sw/nxf_pipelines/telomeric_ends.conf run /lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/sw/nxf_pipelines/telomeric_ends.nf --reads /lustre/scratch123/tol/teams/blaxter/projects/tol-nemotodes/bothBatch/analyses/qualitymetrics/miniBtk-20210123/filteredData/DF5120_filtered.ccs.fasta.gz --outdir testDimiId -profile farm -resume
*/
