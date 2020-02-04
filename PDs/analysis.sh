#!/bin/bash

#The script is editted by Yan-Ru Ju, 2020-1-4
root_dir="/data/usrhome/LabSPLin/splin01/Human_exosomal_miRNA"
seq_dir="$root_dir/src/Human_exosomal_miRNA/190523_Yoshioka_Full/KW_cleanup/processedFiles";
genome_dir="$root_dir/src/hg38";
genomeFa="$genome_dir/hg38.fa";
bt_pref="$root_dir/src/hg38/bt-index/hg38";
miR_db="$root_dir/src/mirbaseDB/";
rRNAtRNA_pref="$root_dir/src/rRNAtRNA.bt-index/hg38.rRNAtRNA";
piRNA_bt_pref="$root_dir/src/piR_bt-idx/piR_bt-idx";  
piRNA_mis3_bt_pref="$root_dir/src/piR-mis3_bt-idx/piR-mis3_bt-idx";  #piR-mis3_bt-idx for full known human piRNAs
star_idx="$genome_dir/hg38.STAR_index";
#STAR_param="$root_dir/src/STAR_piRNA_param.txt"; #Old parameters
STAR_param_YR="$root_dir/src/STAR_piRNA_param_YR.txt";
STAR_mis3_param="$root_dir/src/STAR_piRNA-mis3_param.txt";
patients_name="/data/usrhome/LabSPLin/splin01/YR1033/PDs/patients.list";
Input_miRDeep2_remove_path="/data/usrhome/LabSPLin/splin01/Human_exosomal_miRNA/bulk_analysis/small_RNA_analysis"; 

#consolidate files
##seq_pref=$(ls $seq_dir | sed 's/_S.*//g'|uniq|awk -v f=$1 -v t=$2 '{split($1,a,"_")}a[1]>=f && a[1]<=t{print}');
seq_pref=$(cat $patients_name); ##YR

#run mirdeep mapper.pl
for i in ${seq_pref}; do
  #create analysis directory

  if [ -d "$i" ]; then 
    echo "Removing old folder..."; 
    rm -fR "$i"; 
  fi
  mkdir $i;


  cd $i;
  #concatinate multiple fastq
  echo "Concatinate $i ...";

  if [ -f "$i.fastq" ]; then rm "$i.fastq"; fi
  for f in $(ls $seq_dir|grep $i); do
    echo "Unzipping $f ...";
    gunzip -c $seq_dir/$f >> "$i.fastq";
  done

  #calculate length distribution
  awk 'NR%4==2{s[length($1)]+=1}END{for(i in s){print i"\t"s[i]}}' "$i.fastq" | sort -V > "$i.lenDist";


  #remove rRNAtRNA
  wait;
  echo "Check rRNAtRNA ...";
  bowtie --norc -S -v 0 --al rRNAtRNA.fastq --un rRNAtRNA_removed.fastq $rRNAtRNA_pref tmp.fastq rRNAtRNA.sam > rRNAtRNA.log 2>&1;

  #mirdeep2 mapper
  wait;
  echo "Prepare miRDeep2 input (mapper.pl) for $i ...";
  mapper.pl rRNAtRNA_removed.fastq -n -e -h -j -l 18 -m -p $bt_pref -s reads_collapsed.fa -t reads_vs_refdb.arf -v 2>mapper.log; # -u -k AACTGTAGGCACCATCAAT 

  # mirdeep2
  wait;
  echo "miRDeep2 on $i.fastq ..."; 
  miRDeep2.pl reads_collapsed.fa $genomeFa reads_vs_refdb.arf $miR_db/hsa/mature_hsa.fa $miR_db/hsa/mature_non-hsa.fa $miR_db/hsa/hairpin_hsa.fa -t hsa 2>miRDeep2.log;

  # precursor prediction
  wait;
  echo "Process precursor for $i ...";
  cat mirdeep_runs/run_*/output.mrd | awk '$1~/^>/{nm=$1}$1=="pri_seq"&&nm!="NULL"{print nm"\n"toupper($2);}' | tr "U" "T" > miRDeepPredictedPrecursors.fasta;
  bowtie-build miRDeepPredictedPrecursors.fasta miRDeepPredictedPrecursors;

  # filter out predicted miRNA
  wait;
  echo "Filter with miRNA precursor ...";
  bowtie -norc -S -v 0 --al miRBase.fastq --un miRBase_removed.fastq $miR_db/hsa/hairpin_hsa rRNAtRNA_removed.fastq miRBase.sam > miRBase.log 2>&1;

  wait;
  echo "Mapping with predicted miRNA precursor";
  bowtie -norc -S -v 0 --al miRDeep2.fastq --un miRDeep2_removed.fastq miRDeepPredictedPrecursors miRBase_removed.fastq miRDeep2.sam > miRDeep2.log 2>&1;

  # pick reads within piRNA length range 
  wait;
  echo "Pick reads within piRNA length for analysis";
  awk -v m=26 -v M=34 'BEGIN{OFS="\n"}{nm[NR%4]=$0}NR%4==0{l=length(nm[2]);outLine=nm[1]"\n"nm[2]"\n"nm[3]"\n"nm[0]; if(l>=m && l<=M ){print outLine > "piRNA_candidates.fastq"} else {print outLine > "non-piRNA_candidates.fastq"} }' miRDeep2_removed.fastq;

  # map to known human piRNAs
  wait;
  echo "Mapping for known human piRNAs";
  bowtie -norc -S -v 0 --al piRNA.fastq --un piRNA_removed.fastq $piRNA_bt_pref piRNA_candidates.fastq piRNA.sam > piRNA.log 2>&1;
  bowtie -norc -S -v 0 --al piRNA-mis3.fastq --un piRNA-mis3_removed.fastq $piRNA_mis3_bt_pref piRNA_candidates.fastq piRNA-mis3.sam > piRNA-mis3.log 2>&1;


mkdir $i;
cd $i;

  # map to genome 
  wait;
  echo "Mapping to genome"
  if [ -d "STAR.genome" ]; then rm -fR "STAR.genome"; fi
  if [ -d "unannot.STAR" ]; then rm -fR "unannot.STAR"; fi

  echo "Mapping for post-miRDeep2 reads";
  if [ -d piRNA.STAR ]; then rm -fR piRNA.STAR; fi 
  mkdir piRNA.STAR; 
  ##STAR --genomeDir $star_idx --parametersFiles $STAR_param --outFileNamePrefix ./piRNA.STAR/ --readFilesIn miRDeep2_removed.fastq;

  STAR --genomeDir $star_idx --parametersFiles $STAR_param_YR --outFileNamePrefix ./piRNA.STAR/ --readFilesIn $Input_miRDeep2_remove_path/$i/miRDeep2_removed.fastq; #YR

  wait;


  cd ..
done

