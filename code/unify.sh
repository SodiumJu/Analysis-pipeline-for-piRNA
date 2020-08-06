#!/bin/bash

#Script by Yan-Ru, 2020.5.22

#Intial variables
src_dir="/project/GP1/chialang1220";

sample_id="$src_dir/mingchekuo/Human_fetal_gonad/20130502/sample.name_all.txt";

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

detect_exist $src_dir/Shau-Ping/Human_fetal_gonad/results/$i;

input_main_sam=$src_dir/Shau-Ping/Human_fetal_gonad/process/$i/6_STAR_hg38/Aligned.sortedByCoord.out.sam;
input_rest_sam=$src_dir/Shau-Ping/Human_fetal_gonad/process/$i/6_STAR_hg38/rest.map.STAR.bt.sam;

output_file=$src_dir/Shau-Ping/Human_fetal_gonad/results/$i/$i".uniq.sam";
output_file_bed=$src_dir/Shau-Ping/Human_fetal_gonad/results/$i/$i".uniq.bed12";
output_piRNA_bed=$src_dir/Shau-Ping/Human_fetal_gonad/results/$i/$i".uniq.piRNA.bed12";


awk 'BEGIN{srand();FS=OFS="\t";}
    
       NF<=10{print; next;}
       
       $0~/NH:i:0/ || $6~/I/ || $6~/D/ || $6~/N/ {next;}
       
       {
       
       if($3!~/_/){
            weight=weight+10;
            }
       if(length(hidx[$1])==0){
            hidx[$1]=weight;
            rec[$1]=$0;
            weight=0;
            next;
            }else if(weight>hidx[$1]){
            
            hidx[$1]=weight;
            rec[$1]=$0;
            weight=0;
            next;
             
       }

       }
      
      END{for(i in rec){print rec[i]}}
      ' $input_main_sam > $output_file;

awk 'BEGIN{FS=OFS="\t";}NF<=10{next;}{print;}' $input_rest_sam >> $output_file;

samtools view -Sb $output_file | 
bamToBed -bed12 -i - |
awk '{$4=""; print}' - | sort -V - | uniq -c - |
awk 'BEGIN{OFS="\t"}{ $5=$1; $1=$2; $2=$3; $3=$4; $4="Read_"NR; print; }' - |
tee $output_file_bed |
#output unique reads within piRNA length range
awk -v m=26 -v M=34 'BEGIN{OFS="\t"}
  { split($11,len,","); tLen=0;
  for(l in len){tLen += len[l]} }
  tLen>=m && tLen<=M {print;}' - > $output_piRNA_bed;


done

echo "finished";





