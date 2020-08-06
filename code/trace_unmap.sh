#!/bin/bash
#Script by Yan-Ru, 2020.5.21

#Intial variables

src_dir="/project/GP1/chialang1220";
seq_path="$src_dir/Shau-Ping/Human_fetal_gonad/process/M19-1";

genomeFa="$src_dir/mingchekuo/Reference/UCSU/hg38/hg38.fa";
bt_pref="$src_dir/mingchekuo/Reference/UCSU/hg38/bt-index/hg38";

mkdir $seq_path/7_test_grep_unmap;

#convert files
samtools view -b $seq_path/6_STAR_hg38/unmapped.sam > $seq_path/7_test_grep_unmap/unmapped.bam;
bamToFastq -i $seq_path/7_test_grep_unmap/unmapped.bam -fq $seq_path/7_test_grep_unmap/unmapped.fastq

#bowtie map
cd $seq_path/7_test_grep_unmap;
bowtie --norc -S --no-unal -v 0 --al map.bt.fastq --un removed.bt.fastq $bt_pref $seq_path/7_test_grep_unmap/unmapped.fastq map.bt.sam > test_unmapped.log 2>&1;

echo "done~~";
