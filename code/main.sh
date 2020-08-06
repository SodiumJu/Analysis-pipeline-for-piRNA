#!/bin/bash
#Script by Yan-Ru, 2020.5.13

#Intial variables

src_dir="/project/GP1/chialang1220";
out_dir="$src_dir/Shau-Ping/Human_fetal_gonad/process";
seq_dir="$src_dir/mingchekuo/Human_fetal_gonad/20130502";
sample_id="$src_dir/mingchekuo/Human_fetal_gonad/20130502/sample.name_all.txt";
genomeFa="$src_dir/mingchekuo/Reference/UCSU/hg38/hg38.fa";
bt_pref="$src_dir/mingchekuo/Reference/UCSU/hg38/bt-index/hg38";
miR_db="$src_dir/mingchekuo/Reference/database/mirbase/22";
rRNAtRNA_pref="$src_dir/mingchekuo/Reference/database/GENCODE/human_GRCh38.p12/rRNAtRNA.bt-index/hg38.rRNAtRNA";
piRNA_bt_pref="$src_dir/mingchekuo/Reference/database/piRNAs/piR_bt-idx/piR_bt-idx";  
star_idx="$src_dir/mingchekuo/Reference/UCSU/hg38/hg38.STAR_index";
STAR_param="$src_dir/Shau-Ping/Human_fetal_gonad/process/STAR_piRNA_param.txt";                                                                      

#functions

detect_exist () {
  local File_path=$1;

  if [ -d $File_path ]; then 
    echo "Removing old folder..."; 
    rm -fR $File_path; 
  fi
  mkdir $File_path;
  cd $File_path;
}


seq_pref=$(cat $sample_id);

for i in ${seq_pref}; do
 
  rRNAtRNA_path="$out_dir/$i/1_rRNAtRNA"; #for step 1	
  miRDeep_path="$out_dir/$i/2_miRDeep" #for step 2
  miRNA_path="$out_dir/$i/3_miRNA"; #for step 3
  pre_miRNA_path="$out_dir/$i/4_predicted_miRNA"; #for step 4
  bt_piRdb_path="$out_dir/$i/5_bt_piRdb"; #for step 5
  STAR_hg38_path="$out_dir/$i/6_STAR_hg38"; #for step 6


  detect_exist "$out_dir/$i";
  input_file=$seq_dir/$i"_cutadapt.fastq";

  #Step 0. Calculate length distribution
  echo "Calculate length distribution ...";
  awk 'NR%4==2{s[length($1)]+=1}END{for(i in s){print i"\t"s[i]}}' "$input_file" | sort -V > "$i.lenDist";

  #Step 1. Remove rRNAtRNA
  wait;
  echo "Check rRNAtRNA ...";
  detect_exist "$rRNAtRNA_path";
  bowtie --norc -S -v 0 --al rRNAtRNA.fastq --un rRNAtRNA.removed.fastq $rRNAtRNA_pref $input_file rRNAtRNA.sam > rRNAtRNA.log 2>&1;

  #Step 2. Prepare for mirDeep2
  #mirdeep2 mapper
  wait;
  echo "Prepare miRDeep2 input (mapper.pl) for $i ...";
  detect_exist "$miRDeep_path";
  mapper.pl $rRNAtRNA_path/rRNAtRNA.removed.fastq -n -e -h -j -l 18 -m -p $bt_pref -s reads_collapsed.fa -t reads_vs_refdb.arf -v 2>mapper.log; # -u -k AACTGTAGGCACCATCAAT 

  #mirdeep2
  wait;
  echo "miRDeep2 on $i.fastq ..."; 
  miRDeep2.pl reads_collapsed.fa $genomeFa reads_vs_refdb.arf $miR_db/hsa/mature_hsa.fa $miR_db/hsa/mature_non-hsa.fa $miR_db/hsa/hairpin_hsa.fa -t hsa 2>miRDeep2.log;

  # precursor prediction
  wait;
  echo "Process precursor for $i ...";
  cat mirdeep_runs/run_*/output.mrd | awk '$1~/^>/{nm=$1}$1=="pri_seq"&&nm!="NULL"{print nm"\n"toupper($2);}' | tr "U" "T" > miRDeepPredictedPrecursors.fasta;
  bowtie-build miRDeepPredictedPrecursors.fasta miRDeepPredictedPrecursors;

  #Step 3. Remove miRNA
  # filter out predicted miRNA
  wait;
  echo "Filter with miRNA precursor ...";
  detect_exist "$miRNA_path";
  bowtie --norc -S -v 0 --al miRBase.fastq --un miRBase.removed.fastq $miR_db/hsa/hairpin_hsa $rRNAtRNA_path/rRNAtRNA.removed.fastq miRBase.sam > miRBase.log 2>&1;


  #Step 4. Remove predicted-miRNA
  wait;
  echo "Mapping with predicted miRNA precursor";
  detect_exist "$pre_miRNA_path";
  bowtie --norc -S -v 0 --al miRDeep2.fastq --un miRDeep2.removed.fastq $miRDeep_path/miRDeepPredictedPrecursors $miRNA_path/miRBase.removed.fastq miRDeep2.sam > miRDeep2.log 2>&1;


  #Step 5. Map to human piRNA database (Bowtie)

  # pick reads within piRNA length range 
  wait;
  echo "Pick reads within piRNA length for analysis";
  detect_exist "$bt_piRdb_path";
  awk -v m=26 -v M=34 'BEGIN{OFS="\n"}{nm[NR%4]=$0}NR%4==0{l=length(nm[2]);outLine=nm[1]"\n"nm[2]"\n"nm[3]"\n"nm[0]; if(l>=m && l<=M ){print outLine > "bt_piRdb.fastq"} else {print outLine > "non-bt_piRdb.fastq"} }' $pre_miRNA_path/miRDeep2.removed.fastq;

  # map to known human piRNAs
  wait;
  echo "Mapping for known human piRNAs";
  bowtie --norc -S --no-unal -v 0 --al piRNA.piRdb.fastq --un piRNA.removed.piRdb.fastq $piRNA_bt_pref bt_piRdb.fastq piRNA.piRdb.sam > piRNA.log 2>&1;



  #Step 6. Map to human genome (STAR)
  wait;
  echo "Mapping to genome";
  detect_exist "$STAR_hg38_path";
  STAR --genomeDir $star_idx --parametersFiles $STAR_param --outFileNamePrefix ./ --readFilesIn $pre_miRNA_path/miRDeep2.removed.fastq;
  samtools view -h Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.sam; 
  awk 'BEGIN{FS=OFS="\t";} NF<=10{print; next;} $0~/NH:i:0/||$6~/N/ {print; next;}' Aligned.sortedByCoord.out.sam > STAR.unmap.sam;

  # convert files
  samtools view -b STAR.unmap.sam > STAR.unmap.temp.bam;
  bamToFastq -i STAR.unmap.temp.bam -fq STAR.unmap.temp.fastq;

  # using bowtie map to genome again for those reads that STAR can not handle.
  bowtie --norc -S --no-unal -v 0 --un unmap.STAR.bt.fastq $bt_pref STAR.unmap.temp.fastq rest.map.STAR.bt.sam > unmapped.bt.log 2>&1;
  rm *.temp.*;


done
echo "Finished!!";



