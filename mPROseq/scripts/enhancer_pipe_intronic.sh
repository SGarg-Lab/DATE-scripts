#!/bin/sh
## arg 1 - mm10RefFlat_glist
## arg 2 - path to enhancer scripts
## arg 3 - filtered gene list
## arg 4 - dREG HD output (1+)



a=0
merge="python merge_enhancers.py"
for arg in "$@"
do
    let "a += 1"
    if [ "$a" -gt 3 ]
    then
        python $2/enhancer_filter_intronic.py $arg $1 > ${arg%.*}_enhancers_unsorted.txt
        sort -k 1,1 -k 2,2n -k 3,3n ${arg%.*}_enhancers_unsorted.txt > ${arg%.*}_enhancers.txt
        merge="$merge ${arg%.*}_enhancers.txt"
    fi
done

if [ "$a" -lt 5 ]
then
    python $2/arrange.py ${4%.*}_enhancers.txt > ${4%.*}_arranged.txt

    awk 1 $3 ${4%.*}_arranged.txt > ${4%.*}_enhgene_unsorted_tmp.txt

    awk 'FNR>1' ${4%.*}_enhgene_unsorted_tmp.txt > ${4%.*}_enhgene_unsorted.txt

    sort -k 3,3 -k 5,5n -k 6,6n ${4%.*}_enhgene_unsorted.txt > ${4%.*}_enhgene_unnum.txt

    python $2/add_nums.py ${4%.*}_enhgene_unnum.txt > ${4%.*}_kbro_input.txt

    python $2/stats_filter.py ${4%.*}

    mkdir ${4%.*}_enhancerlists
    mv ${4%.*}_dropped.txt ${4%.*}_enhancerlists
    mv ${4%.*}_enhancers.txt ${4%.*}_enhancerlists
    mv ${4%.*}_kbro_input.txt ${4%.*}_enhancerlists
    mv ${4%.*}_dist-size.png ${4%.*}_enhancerlists
    mv $4 ${4%.*}_enhancerlists

    rm ${4%.*}_enhancers_unsorted.txt
    rm ${4%.*}_enhgene_unsorted_tmp.txt
    rm ${4%.*}_arranged.txt
    rm ${4%.*}_enhgene_unnum.txt
    rm ${4%.*}_enhgene_unsorted.txt
else
    $merge

    python $2/arrange.py enhancers_merged.bed > enhancers_arranged.txt

    awk 1 $1 enhancers_arranged.txt > enhancers_enhgene_unsorted.txt

    sort -k 3,3 -k 5,5n -k 6,6n enhancers_enhgene_unsorted.txt > enhancers_enhgene_unnum.txt

    python $2/add_nums.py enhancers_enhgene_unnum.txt > enhancers_kbro_input.txt

    python $2/stats_filter.py enhancers

    rm enhancers_unsorted.txt
    rm enhancers_arranged.txt
    rm enhancers_enhgene_unnum.txt
    rm enhancers_enhgene_unsorted.txt
fi
