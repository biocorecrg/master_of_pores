#!/bin/bash

if [ x"$1" == x ]; then 
        echo "please specify a bam file"
        exit 1 
fi

if [ x"$2" == x ]; then 
        echo "please specify a gtf file with annotations! The tag gene_id must be present"
        exit 1 
fi

if [ x"$3" == x ]; then 
        echo "please specify an output file"
        exit 1 
fi

intersectBed -b $1 -a $2 -bed -wb -s |\
awk -F "\t" '{if ($3=="gene") {
	split($9,a,";"); for (i in a) 
		{ if (a[i] ~ "gene_id") 
			{ gsub(/gene_id/, "", a[i]); gsub(/"/, "", a[i]); print $13"\t"a[i] 
			}
		}
	}
}' > $3
