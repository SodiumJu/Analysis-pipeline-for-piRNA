#!/bin/bash

#PBS -P MST108173
#PBS -W group_list=MST108173
#PBS -N R403_piRNA_test
#PBS -l select=1:ncpus=40
#PBS -l place=pack
#PBS -q ngs192G
#PBS -V
#PBS -o /project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code/logs/ 
#PBS -e /project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code/logs/

local_path="/project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code";

source $HOME/.bashrc;
#conda activate ncRNA;

#$local_path/main.sh;
#wait;
#$local_path/unique.mapping.sh;

conda deactivate;

conda activate piRNA;

$local_path/runPiRWeblogo.sh
#$local_path/trace_unmap.sh;

echo "PBS success";

