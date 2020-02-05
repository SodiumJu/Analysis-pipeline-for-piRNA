#!/bin/bash

oFile="TE_summary.csv";

#build header line
#items=("sample_id" "Total_reads" "Genic" "Repeat" "Genic&Repeat" "Intergenic");
items=("Sample_ID" "SINE?" "Satellite" "rRNA" "snRNA" "DNA?" "scRNA" "tRNA" "DNA" "RNA" "RC?" "SINE" "Unknown" "RC" "LTR?" "Low_complexity" "LTR" "Retroposon" "Simple_repeat" "LINE" "srpRNA");
#printf '%s\n' ${items[@]} | paste -sd ',' >> "$oFile";
echo "${items[@]}" > "$oFile";

# f and t are internal control indicate file range for output
for i in $(ls|grep "repeatmasker.Class" |sed 's/.repeatmasker.Class//g'); do
  echo "Processing $i ...";
  SINE_p=`cat $i.repeatmasker.Class|grep 'SINE?'|sed 's/SINE? //g'`;
  Satellite=`cat $i.repeatmasker.Class|grep 'Satellite'|sed 's/Satellite //g'`;
  rRNA=`cat $i.repeatmasker.Class|grep 'rRNA'|sed 's/rRNA //g'`;
  snRNA=`cat $i.repeatmasker.Class|grep 'snRNA'|sed 's/snRNA //g'`;
  DNA_p=`cat $i.repeatmasker.Class|grep 'DNA?'|sed 's/DNA? //g'`;
  scRNA=`cat $i.repeatmasker.Class|grep 'scRNA'|sed 's/scRNA //g'`;
  tRNA=`cat $i.repeatmasker.Class|grep 'tRNA'|sed 's/tRNA //g'`;
  DNA=`cat $i.repeatmasker.Class|grep '^DNA '|sed 's/DNA //g'`;
  RNA=`cat $i.repeatmasker.Class|grep '^RNA '|sed 's/RNA //g'`;
  RC_p=`cat $i.repeatmasker.Class|grep 'RC?'|sed 's/RC? //g'`;
  SINE=`cat $i.repeatmasker.Class|grep 'SINE '|sed 's/SINE //g'`;
  Unknown=`cat $i.repeatmasker.Class|grep 'Unknown'|sed 's/Unknown //g'`;
  RC=`cat $i.repeatmasker.Class|grep 'RC '|sed 's/RC //g'`;
  LTR_p=`cat $i.repeatmasker.Class|grep 'LTR?'|sed 's/LTR? //g'`;
  Low_complexity=`cat $i.repeatmasker.Class|grep 'Low_complexity'|sed 's/Low_complexity //g'`;
  LTR=`cat $i.repeatmasker.Class|grep 'LTR '|sed 's/LTR //g'`;
  Retroposon=`cat $i.repeatmasker.Class|grep 'Retroposon '|sed 's/Retroposon //g'`;
  Simple_repeat=`cat $i.repeatmasker.Class|grep 'Simple_repeat'|sed 's/Simple_repeat //g'`;
  LINE=`cat $i.repeatmasker.Class|grep 'LINE '|sed 's/LINE //g'`;
  srpRNA=`cat $i.repeatmasker.Class|grep 'srpRNA'|sed 's/srpRNA //g'`;
  #Total_reads=`cat $i'_Refproportion'|grep 'Total_reads'|sed 's/Total_reads //g'`;
  #Genic=`cat $i'_Refproportion'|grep '^Genic '|sed 's/Genic //g'`;
  #Repeat=`cat $i'_Refproportion'|grep '^Repeat'|sed 's/Repeat //g'`;
  #Genic_and_Repeat=`cat $i'_Refproportion'|grep 'Genic&Repeat'|sed 's/Genic&Repeat //g'`;
  #Intergenic=`cat $i'_Refproportion'|grep 'Intergenic'|sed 's/Intergenic //g'`;
  echo "Input $i done";
  samp=("$i" "$SINE_p" "$Satellite" "$rRNA" "$snRNA" "$DNA_p" "$scRNA" "$tRNA" "$DNA" "$RNA" "$RC_p" "$SINE" "$Unknown" "$RC" "$LTR_p" "$Low_complexity" "$LTR" "$Retroposon" "$Simple_repeat" "$LINE" "$srpRNA");
  echo ${samp[@]} >> "$oFile";
done

#substitute separater to tab
sed -i 's/ /,/g' "$oFile";
