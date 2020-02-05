#!/bin/bash

#Program:
#Calculating piRNAs originated from TE belonging to which classes
#History:
#2018/11/07 HHL
#Modified by YR in 2019/11/4
#Because 5th column is the number of piRNA counts and 17th column is the class name of TEs
 
#inputdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/bed/
inputdir=/data/usrhome/LabSPLin/splin01/Human_exosomal_miRNA/bulk_analysis/piRNA_analysis/uniqPirBed;
refdir=$HOME/genomes/UCSC/hg38/RepeatMasker;
#refdir=$HOME/YR1033/PDs/TEclass_proportion;
for input in $(ls $inputdir|grep '.bed12'|sed 's/.uniq.piR.bed12//g')
do

echo Start aligning $input to repeat region and counting reads;
intersectBed -s -wa -wb -a $inputdir/$input.uniq.piR.bed12 -b $refdir/hg38_RepeatMasker_annotated.bed > repeatmasker.temp;

echo Start classifing TEs;
awk '{split ($17,a,"|");print a[2],"\t",$5}' repeatmasker.temp|awk '{a[$1]=a[$1]?a[$1]+$2:$2;}END{for (i in a)print i, a[i];}'>$input.repeatmasker.Class; #YR sum of each catagories
#echo Start classifing TEs;

#cat repeatmasker.temp|awk '{split ($16,a,"|");print a[2] $5;}'|sort -V|uniq -c|sort -r -n > $input.repeatmasker.Class;

echo $input is finished;
rm repeatmasker.temp
done

