#!/bin/bash

#Yan-Ru 2020.5.28
 
doWeblogoPNG () {
  local inFile=$1;
  local outName=$2;

  if [[ "$inFile" =~ ".bed12" ]]; then
    echo converting condensed bed12 to fasta;
    cat $inFile |
      # expand to actual copy of reads
      awk 'BEGIN{OFS="\t"}{o=$4;for(i=0;i<$5;i++){$4=o"_"i; print}}' - |
      # compute fasta
      fastaFromBed -name+ -s -split -fi $fa_path -bed - -fo $outName".tmp.fasta";
  else
      cp $inFile $outName".tmp.fasta";
  fi
  
  cat $outName".tmp.fasta" | 
    # format into RNA (replace T with U)
    tr 'tT' 'U' |
    awk -v readLen=20 'NR%2==1{print}NR%2==0{ print substr($1,1,readLen) }' - > $outName".weblogo.tmp";

  for unit in "bits" "probability"; do
    echo $outName"."$unit".png"; 
    weblogo -f $outName".weblogo.tmp" -o $outName"."$unit".png" -F png_print -A rna -U $unit --fineprint seqlogo -c classic;
    wait;
  done
  rm $outName".tmp.fasta" $outName".weblogo.tmp";
}

find10bpOverlap () {
  local input_file=$1
  local Name=$2;
  local outpath=$3;
  echo "Finding 10bp overlap for $pfx$inBed12$sfx ..." ;
  echo "Outputting $inBed12.10bpOverlap.bed12 ...";
  awk -v out=$outpath '$6=="+"{print > out"/p.bed12"}$6=="-"{print > out"/n.bed12"}' $input_file;
  intersectBed -wo -S -a $outpath/p.bed12 -b $outpath/n.bed12 | 
    awk '$3>$15 && $15>$2 && $2>$14 && $15-$2==9{
           for(i=1;i<NF;i++) {
             printf "%s", $i; 
             if(i%12>0){printf "\t"} else{printf "\n"}
         } }' |
    sort -V | uniq > $outpath/$Name".10bpOverlap.bed12"; 
  #rm $outpath/p.bed12 $outpath/n.bed12;
}

findTE () {
  local input_file=$1;
  local Name=$2;
  local outpath=$3;
  echo "Finding TE for $pfx$inBed12$sfx ..." ;
  echo "Outputting $inBed12'.TE.bed12' ...";
  intersectBed -split -u -wa -a $input_file -b $TEdir -f 0.7 > $outpath/$Name".TE.bed12";
}

# Parameters 
srcDir="/project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad"
subDir="$srcDir/results";
subDir_process="$srcDir/process";
fInSfx=( ".uniq.piRNA.bed12" "piRNA.piRdb.sam");

TEdir="/project/GP1/chialang1220/mingchekuo/Reference/UCSU/hg38/Repeatmasker/hg38_RepeatMasker.bed";
outDir="$srcDir/results";
sampName="/project/GP1/chialang1220/mingchekuo/Human_fetal_gonad/20130502/sample.name_all.txt";
fa_path="/project/GP1/chialang1220/mingchekuo/Reference/UCSU/hg38/hg38.fa";

# main

#process bed12 files
#for length selection putative piRNA

samples=$(cat $sampName);
for fName in ${samples} ; do
  
  #create sub-directories

  for outSubDir in "full" "10bp" "TEassoc" "TEassoc_10bp" "piRdb" "piRdb_POS1"; do
	 
      if [ -d $outDir/$fName/$outSubDir ]; then 
		      
			  echo "Removing old folder..."; 
			  rm -fR $outDir/$fName/$outSubDir;
			  mkdir $outDir/$fName/$outSubDir; 
	  fi
  done
  
  echo "Identifying 10bpOverlap piRNAs for $fName${fInSfx[0]} ...";
  find10bpOverlap "$subDir/$fName/$fName${fInSfx[0]}" "$fName" $outDir/$fName/10bp;

  echo "Identifying TE associated piRNAs for $fName${fInSfx[0]} ...";
  findTE "$subDir/$fName/$fName${fInSfx[0]}" "$fName" $outDir/$fName/"TEassoc";

  echo "Identifying TE associated 10bpOverlap piRNAs for $fName${fInSfx[0]} ...";
  findTE $outDir/$fName/10bp/$fName".10bpOverlap.bed12" "$fName" $outDir/$fName/"TEassoc_10bp";

  doWeblogoPNG $subDir/$fName/$fName".uniq.piRNA.bed12" $outDir/$fName/full/"full_piR."$fName".weblogo";
  doWeblogoPNG $subDir/$fName/10bp/$fName".10bpOverlap.bed12" $outDir/$fName/10bp/"10bp_piR."$fName".weblogo";
  doWeblogoPNG $subDir/$fName/TEassoc/$fName".TE.bed12" $outDir/$fName/TEassoc/"TEassoc_piR."$fName".weblogo";
  doWeblogoPNG $subDir/$fName/TEassoc_10bp/$fName".TE.bed12" $outDir/$fName/TEassoc_10bp/"TEassoc_10bp_piR."$fName".weblogo";

done

#concatenate fasta files
#for human piRNA database mapped piRNA
<< 'MULTILINE-COMMENT'
for fName in ${samples} ; do
  echo "Processing $fName ...";
  cat $subDir_process/$fName/5_bt_piRdb/${fInSfx[1]} | samtools view -Sb - | samtools fasta -G 4 - > $subDir/$fName/tmp.fa;
  doWeblogoPNG $subDir/$fName/tmp.fa $outDir/piRdb/"hg38_piR."$fName".weblogo";
  
  wait;
  awk 'NF<5 || (NF>5 && $4==1)' $subDir_process/$fName/5_bt_piRdb/${fInSfx[1]} | samtools view -Sb - | samtools fasta -G 4 - > $subDir/$fName/tmp.fa;
  doWeblogoPNG $subDir/$fName/tmp.fa $outDir/piRdb_POS1/"hg38_piR_POS1."$fName".weblogo";
  
  wait;
  rm $subDir/$fName/tmp.fa;
done
MULTILINE-COMMENT


