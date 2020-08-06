#!/bin/bash

oFile="summary.log";

# label sample types
sampTypeFile="sample_type.txt";
declare -A sampType;
while read tmp1 tmp2; do
  sampType+=( [$tmp1]=$tmp2 );
done < $sampTypeFile; 

#build header line
items=( "sample_id" "type" "input_reads" "mappable_reads" "rRNAtRNA" "miRNA" "pred_miRNA" 
	"hsa_piRNA(#POS1)_perfMatch_hg38" "hsa_piRNA(#POS1)_full_piRDB"
	"mappable_unannotated" "mappable_within_piR_length" "unmappable");
echo "${items[@]}" > "$oFile";

# f and t are internal control indicate file range for output
for i in $(ls | grep "_" |awk -v f=$1 -v t=$2 '$1~/_/{split($1,nm,"_")}nm[1]>=f&&nm[1]<=t{print}'); do
  echo "Processing $i ..."; 
  tot_cnt=`cat "$i/rRNAtRNA.log" | awk '$0~/reads processed/{print $NF}'` ;
  rRNAtRNA_cnt=`cat "$i/rRNAtRNA.log" | awk '$0~/reported alignment/{print $(NF-1)}'` ;
  miRBase_cnt=`cat "$i/miRBase.log"  | awk '$0~/reported alignment/{print $(NF-1)}'` ;
  miRDeep2_cnt=`cat "$i/miRDeep2.log" | awk '$0~/reported alignment/{print $(NF-1)}'` ;
  hsa_piR_cnt=`cat "$i/piRNA.log"    | awk '$0~/reported alignment/{print $(NF-1)}'` ;
  piR_st1_cnt=`cat "$i/piRNA.sam" | awk '$4==1{c+=1}END{print c}'` ;
  hsa_piR_mis3_cnt=`cat "$i/piRNA-mis3.log" | awk '$0~/reported alignment/{print $(NF-1)}'` ;
  piR_mis3_st1_cnt=`cat "$i/piRNA-mis3.sam" | awk '$4==1{c+=1}END{print c}'` ;
  mapped_piR_size_cnt=`samtools view -F "0x4" -F "0x100" $i/piRNA.STAR/Aligned.sortedByCoord.out.bam | awk -v m=26 -v M=34 '{l=length($10)}l>=m&&l<=M{s++}END{print s}'`;
  unannot_cnt=`cat "$i/piRNA.STAR/Log.final.out" |
	       awk -v piR=$mapped_piR_size_cnt '
		    $0~/UNIQUE READS/||$0~/MULTI-MAPPING READS/{mappable=1;next}
		    mappable==1{mapped+=$NF; mappable=0;}
		    $0~/UNMAPPED READS/{unmap=1;next;}
		    unmap==1&&$NF!~/%/{unmapped+=$NF;}
		    END{printf "%d %d %d\n", mapped, piR, unmapped;}'` ;
  mappable=`echo $unannot_cnt | awk -v t=$tot_cnt '{print t-$3}'`;
  samp=("$i" "${sampType[$i]}" "$tot_cnt" "$mappable" "$rRNAtRNA_cnt" "$miRBase_cnt" "$miRDeep2_cnt" 
	"$hsa_piR_cnt($piR_st1_cnt)" "$hsa_piR_mis3_cnt($piR_mis3_st1_cnt)" "$unannot_cnt");
  echo ${samp[@]} >> $oFile; 
done

#substitute separater to tab
sed -i 's/ /\t/g' "$oFile";


