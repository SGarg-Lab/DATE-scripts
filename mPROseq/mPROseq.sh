#!/bin/bash

# example script for processing mPROseq reads


module load fastxtoolkit/0.0.13 bowtie2/2.3.5.1 samtools/1.10
module load python/2.7.13 deeptools/3.0.1
module load bedtools/2.29.2 ucsc-tools/20170321

# directory contains fastq files and text file listing sample index for each library
cd /path/to/directory/

# names of fastq files to analyze
List="
mPROseq_poolA
mPROseq_poolB
"

mm10=path/to/bowtie2/index/
dREG_model=path/to/asvm.gdm.6.6M.20170828.rdata #downloaded from https://cbsuftp.tc.cornell.edu/danko/hub/dreg.models/asvm.gdm.6.6M.20170828.rdata
dREG_HD_model=path/to/dREG_HD.model.rdata # downloaded from https://github.com/Danko-Lab/dREG.HD/blob/master/dREG.HD/inst/extdata/dREG_HD.model.rdata
chromSizesmm10=path/to/mm10.chrom.sizes # downloaded from https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
glist1=path/to/refFlat.txt # downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
glist2=path/to/glist_mm10_NoAltTranscripts_ReAnnotatedTSS_UniqueName_Jay.txt
script=path/to/scripts

mkdir split_fasta

for x in ${List}; do
    echo converting fastq to fasta file from $x:
    fastq_to_fasta -n -v -i $x.fastq -o $x.fasta -Q33

    echo clipping adapter sequence from $x:
    fastx_clipper -v -n -a ATCTCGTATGCCGTCTTCTGCTTG -l 27 -i $x.fasta -o ${x}_adapclip.fasta -Q33

    echo clipping adapter sequence from $x:
    fastx_clipper -v -n -a GGGGGGGGGGGG -l 27 -i ${x}_adapclip.fasta -o ${x}_adapclip2.fasta -Q33

    echo collapsing reads using UMIs from $x:
    fastx_collapser -v -i ${x}_adapclip2.fasta -o ${x}_collapsed.fasta -Q33

    echo trimming UMI from $x:
    fastx_trimmer -f 7 -v -i ${x}_collapsed.fasta -o ${x}_trimmed.fasta -Q33

    echo barcode-splitting the reads from $x using sample index:
    cat ${x}_trimmed.fasta | fastx_barcode_splitter.pl --bcfile Index_$x.txt --bol --prefix split_fasta/ --mismatches 3 --partial 2 --suffix ".fasta"
done

cd split_fasta
mkdir bed bedgraph bedgraph_RPMnorm bigWig
for x in *.fasta; do
  label=$(echo $x | awk '{gsub(".fasta", "")}1')

  echo trimming the sample index and additional C in front of sample index for $label:
  fastx_trimmer -Q33 -f 10 -v -i ${label}.fasta -o ${label}_trimmed.fasta

  echo finding reverse complement for $in:
  fastx_reverse_complement -Q33 -i ${label}_trimmed.fasta -o ${label}_RC.fasta

  echo aligning reads to mm10 for $in:
  bowtie2 --sensitive-local -p 7 -x $mm10 -U ${label}_RC.fasta -s '-' | samtools view -S -b '-' > ${label}.bam
  samtools sort ${label}.bam -o ${label}_sorted.bam
  samtools index ${label}_sorted.bam

  echo converting bam to bed file and sorting it:
	bedtools bamtobed -i ${label}.bam | sort -k1,1 -k2,2n > bed/${label}_sorted.bed

  echo generating non-normalized bedgraph files of ${label}:
	awk '$6 == "-"' bed/${label}_sorted.bed | bedtools genomecov -i stdin -3 -bg -g ${chromSizesmm10} > bedgraph/${label}_m.bedgraph
	awk '{$4=$4*-1; print}' bedgraph/${label}_m.bedgraph > bedgraph/${label}_mn.bedgraph
  awk '$6 == "+"' bed/${label}_sorted.bed | bedtools genomecov -i stdin -3 -bg -g ${chromSizesmm10} > bedgraph/${label}_pl.bedgraph

  echo generating RPM normalized bedgraphs files of ${label}:
	a=$(awk '{ sum += $4 } END { print sum }' bedgraph/${label}_pl.bedgraph )
	b=$(awk '{ sum += $4 } END { print sum }' bedgraph/${label}_mn.bedgraph )
	d=$(($b*-1))
	c=$(expr $a + $d)

	echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' bedgraph/${label}_pl.bedgraph > bedgraph_RPMnorm/${label}_RPMnorm_pl.bedgraph
	echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' bedgraph/${label}_mn.bedgraph > bedgraph_RPMnorm/${label}_RPMnorm_mn.bedgraph

  echo making bigwig from normalized bedgraph files of ${label}:
  bedGraphToBigWig bedgraph_RPMnorm/${label}_RPMnorm_pl.bedgraph ${chromSizesmm10} bigWig/${label}_RPMnorm_pl.bigWig
  bedGraphToBigWig bedgraph_RPMnorm/${label}_RPMnorm_mn.bedgraph ${chromSizesmm10} bigWig/${label}_RPMnorm_mn.bigWig

  echo running dREG:
  $script/run_dREG_HD.bsh bigWig/${label}_RPMnorm_pl.bigWig bigWig/${label}_RPMnorm_mn.bigWig $label $dREG_model 14

  echo running dREG HD:
  gunzip $label.dREG.peak.score.bed.gz > $label.dREG.peak.score.bed
  $script/run_DREG_HD.bsh $label.dREG.peak.score.bed bigWig/${label}_RPMnorm_pl.bigWig bigWig/${label}_RPMnorm_mn.bigWig $dREG_HD_model 16
  mv $label.dREG.peak.score.bed_dREG_HD_stringent.bed ${label}_dREG_HD.bed

  echo building enhancer list:
  $script/enhancer_pipe_intronic.sh $glist1 $glist2 ${label}_dREG_HD.bed

done
