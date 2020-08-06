#!/bin/bash         

bt_pref="$src_dir/mingchekuo/Reference/UCSU/hg38/bt-index/hg38";
 M19.unmap.STAR.bt.fastq
pref="/project/GP1/chialang1220/Shau-Ping/Human_fetal_gonad/code/test";
M19-1_piRNA_candidate.fastq
bowtie --norc -S --no-unal -v 0 --un unmap.STAR.bt.fastq $bt_pref STAR.unmap.temp.fastq rest.map.STAR.bt.sam > unmapped.bt.log 2>&1;
