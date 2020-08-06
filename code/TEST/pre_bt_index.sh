#!/bin/bash         

$1="M19-1_piRNA_candidate";

#fastq to fasta
#sed -n '1~4s/^@/>/p;2~4p' "$1".fastq > "$1".fa;

bowtie-build "$1".fa "$1";
