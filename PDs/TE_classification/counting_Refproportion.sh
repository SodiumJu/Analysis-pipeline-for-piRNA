#!/bin/bash

#Program:
#	Calculating the proportion of piRNA originated from different genome region

#History
#2019/10/21 YR1033
#	Modified for human plasma small RNAs


inputdir=/data/usrhome/LabSPLin/splin01/Human_exosomal_miRNA/bulk_analysis/piRNA_analysis/uniqPirBed
refdir=/data/usrhome/LabSPLin/splin01/genomes/UCSC/hg38 #YR

for input in $(ls $inputdir/|grep '.bed12'|sed 's/.uniq.piR.bed12//g') #YR
do

file_name=$input'_Refproportion'; #YR
echo Input $input small RNA-seq data;
echo "" > "$file_name"; #YR
echo --$input-- >> "$file_name"; #YR;

#Counting total reads
total=`awk '{x+=$5;}END{print x}' $inputdir/$input.uniq.piR.bed12`; #YR
echo Total_reads $total >> "$file_name"; #YR
wait;
#total=`cat $inputdir/$input.bed12|wc -l`;

#Counting reads originated from genic region

echo Start aligning $input to genic region and counting reads;
intersectBed -s -u -wa -a $inputdir/$input.uniq.piR.bed12 -b $refdir/RefSeq/hg38_RefSeq_annotated.bed > refseq.temp; #YR
genic_fusion=`awk '{x+=$5;}END{print x}' refseq.temp`; #YR
wait;

#Counting reads originated from repeat region
echo Start aligning $input to repeat region and counting reads;
intersectBed -s -u -wa -a $inputdir/$input.uniq.piR.bed12 -b $refdir/RepeatMasker/hg38_RepeatMasker_annotated.bed > repeatmasker.temp; #YR
repeat_fusion=`awk '{x+=$5;}END{print x}' repeatmasker.temp`; #YR
wait;
#repeat=`cat repeatmasker.temp|wc -l`;

#Counting reads originated from genic&repeat fusion region 
echo Start aligning $input to genic"&"repeat fusion region and counting reads;
intersectBed -s -u -wa -a refseq.temp -b $refdir/RepeatMasker/hg38_RepeatMasker_annotated.bed > fusion.temp; #YR
fusion=`awk '{x+=$5;}END{print x}' fusion.temp`; #YR
wait;
#fusion=`cat fusion.temp|wc -l`;

#Counting reads #YR
genic=$(($genic_fusion-$fusion));
repeat=$(($repeat_fusion-$fusion));
intergenic=$(($total-$genic-$repeat-$fusion));

#Recording reads #YR
echo Start counting reads from intergenic region;
echo Genic $genic >> "$file_name";
echo Repeat $repeat >> "$file_name";
echo Genic"&"Repeat $fusion >> "$file_name";
echo Intergenic $intergenic >> "$file_name";
wait;


#echo Start counting reads from intergenic region;
#echo Intergenic $(($total-$genic-$repeat+$fusion)) >> Refproportion;
#echo Reads originated from intergenic region = $(($total-$genic-$repeat+$fusion));

#Finish and remove al the temporary files
echo $input is finished,remove temp;
rm *.temp;
echo "" >> "$file_name"; #YR
wait;
echo $input' is done'; #YR

done




