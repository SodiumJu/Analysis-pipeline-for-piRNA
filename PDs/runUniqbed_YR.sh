#!/bin/bash

#srcDir="$HOME/Human_exosomal_miRNA/bulk_analysis/small_RNA_analysis";
srcDir="$HOME/YR1033/PDs/1_Small_RNA_analysis"; #YR
#fInSfx="piRNA.STAR/Aligned.sortedByCoord.out.bam";
fInSfx="piRNA.STAR/Aligned.sortedByCoord.out.bam"; #YR

outDir="$HOME/YR1033/PDs/2_transform_files/uniqBed"; #YR
fOutSfx="uniq.bed12";

outPiRDir="$HOME/YR1033/PDs/2_transform_files/uniqPirBed"; #YR
fPiROutSfx="uniq.piR.bed12";

patients_name="/data/usrhome/LabSPLin/splin01/YR1033/PDs/patients.list"; #YR

if [ ! -d $outDir ]; then mkdir $outDir; fi
if [ ! -d $outPiRDir ]; then mkdir $outPiRDir; fi

seq_pref=$(cat $patients_name); ##YR

#run mirdeep mapper.pl
for fName in ${seq_pref}; do

#for fName in $(ls $srcDir|grep "_"); do
  echo "Processing $fName";
  # make bed12 format
  samtools view -h $srcDir/$fName/$fInSfx |  
    # retain only first mapping for multimap reads
    awk 'BEGIN{FS=OFS="\t"} 
       NF<=10{print; next;} 
       $0~/NH:i:0/{next;} 
       { for(i=1;i<=NF;i++){ if($i~/HI:i:/){split($i,n,":")} } } 
       length(hidx[$1])==0{rec[$1]=$0; hidx[$1]=idx; hspl[$1]=0; if($6~/N/){hspl[$1]=1}; next;} 
       hidx[$1]>n[3] && $6!~/N/{rec[$1]=$0; hidx[$1]=n[3]; hspl[$1]=0; next;} 
       hidx[$1]>n[3] && hspl[$1]==1{rec[$1]=$0; hidx[$1]=n[3]; hspl[$1]=1; next;} 
       END{for(i in rec){print rec[i]}}
      ' - |
    samtools view -Sb |
    # convert to bed12 format 
    bamToBed -bed12 -i - |
    # "compress" reads mapped to same site
    awk '{$4=""; print}' - | sort -V - | uniq -c - |
    awk 'BEGIN{OFS="\t"}{ $5=$1; $1=$2; $2=$3; $3=$4; $4="Read_"NR; print }' - | 
    #output unique reads
    tee $outDir/$fName.$fOutSfx | 
    #output unique reads within piRNA length range
    awk -v m=26 -v M=34 'BEGIN{OFS="\t"} 
      { split($11,len,","); tLen=0;
        for(l in len){tLen += len[l]} }
      tLen>=m && tLen<=M {print;}' - > $outPiRDir/$fName.$fPiROutSfx; 
done



