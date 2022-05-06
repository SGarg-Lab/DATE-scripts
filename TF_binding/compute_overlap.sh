#!/bin/bash

# example script for computing Jaccard index for TF binding

module load bedtools/2.29.2

cd /path/to/narrowPeak/files
window=500
dates=path/to/dates.bed

mkdir sorted

# convert macs2-generated narrowPeak to bed format
for file in *.narrowPeak; do
   name=$(awk -F_peaks '{print $1}' <<< $file)
   suffix=".bed"
   mv $file $name$suffix
 done

# #sorting by chrm and start
 for filename in *.bed; do
   sort -V -k1,1 -k2,2 -k3,3 $filename > sorted/$filename
 done


cd sorted 
mkdir filesizes DATEs paired_all paired_DATEs

#--output original file lengths per TF
for tf in *.bed; do
  echo $(wc -l $tf) >> filesizes/lengths.txt
done
#--output number of paired all
for tf1 in *.bed; do
  for tf2 in *.bed; do
    bedtools window -a $tf1 -b $tf2 -w $window -u > paired_all/$tf1$tf2
  done
done
#--output dates that overlap TF
for tf in *.bed; do
  bedtools window -a $dates -b $tf -w $window -u > DATEs/$tf
done

#--output dates that overlap both
for tf1 in DATEs/*.bed; do
  for tf2 in DATEs/*.bed; do
    bedtools intersect -a DATEs/$tf1 -b DATEs/$tf2 -f 1 -r -u > paired_DATEs/$tf1$tf2
  done
done

#--output number of overlaps
cd DATEs
for tf1tf2 in *.bed; do
  echo $(wc -l $tf1tf2) >> ../filesizes/datetf.txt
done
cd ../paired_DATEs
for tf1tf2 in *.bed; do
  echo $(wc -l $tf1tf2) >> ../filesizes/datepaired.txt
done
cd ../paired_all
for tf1tf2 in *.bed; do
  echo $(wc -l $tf1tf2) >> ../filesizes/allpaired.txt
done
