#!/bin/bash

#PBS -P MST108173
#PBS -W group_list=MST108173
#PBS -N R403_piRNA_test
#PBS -l select=1:ncpus=20
#PBS -l place=pack
#PBS -q ngs96G
#PBS -V
#PBS -o /project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code/logs/ 
#PBS -e /project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code/logs/

local_path="/project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/process/M19-1/6_STAR_hg38";

source $HOME/.bashrc;
conda activate ncRNA;

./pre_bt_index.sh M19-1_piRNA_candidate;

#samtools view -h $local_path/Aligned.sortedByCoord.out.bam > $local_path/Aligned.sortedByCoord.out.sam;
echo "PBS success";

