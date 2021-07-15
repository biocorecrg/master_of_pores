#!/usr/bin/env bash
#
# This script convert the tombo output in a 6 field bed file t

# specify the name of the tombo input
if [ x"$1" == x ]; then

        echo "please specify the tombo output file" 

        exit 1

fi

if [ x"$2" == x ]; then

        echo "please specify the name of the output file"

        exit 1

fi


grep ">" $1 |  awk '{print $1"\t"$5}' | sed s/">"//g |  sed s/\:/\\t/g |   awk '{num++; OFS="\t"; print $1,$2,$2, $4*100, "PRED_"num,$3}' > $2
