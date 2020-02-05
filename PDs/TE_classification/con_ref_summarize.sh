#!/bin/bash

oFile="summary.csv";

#build header line
items=("sample_id" "Total_reads" "Genic" "Repeat" "Genic&Repeat" "Intergenic");
#printf '%s\n' ${items[@]} | paste -sd ',' >> "$oFile";
echo "${items[@]}" > "$oFile";

# f and t are internal control indicate file range for output
for i in $(ls|grep "[0-9]_Refproportion" |sed 's/_Refproportion//g'); do
  echo "Processing $i ...";
  Total_reads=`cat $i'_Refproportion'|grep 'Total_reads'|sed 's/Total_reads //g'`;
  Genic=`cat $i'_Refproportion'|grep '^Genic '|sed 's/Genic //g'`;
  Repeat=`cat $i'_Refproportion'|grep '^Repeat'|sed 's/Repeat //g'`;
  Genic_and_Repeat=`cat $i'_Refproportion'|grep 'Genic&Repeat'|sed 's/Genic&Repeat //g'`;
  Intergenic=`cat $i'_Refproportion'|grep 'Intergenic'|sed 's/Intergenic //g'`;
  echo "Input $i done";
  samp=("$i" "$Total_reads" "$Genic" "$Repeat" "$Genic_and_Repeat" "$Intergenic");
  echo ${samp[@]} >> "$oFile";
done

#substitute separater to tab
sed -i 's/ /,/g' "$oFile";
