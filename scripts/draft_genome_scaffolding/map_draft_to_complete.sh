#!/bin/bash
references="./references/*.fna"
for assembly in `ls -1 *.fna`;  do
	dataset="${assembly/.fna/}"
	for reference in $references; do
		outname=`echo $reference | sed 's/.fna//g;s/^.*\///g'` 
		python2.7 ~/Software/CONTIGuator_v2.7/CONTIGuator.py -r $reference -c $assembly -M -V -f ${dataset}_$outname -t 6 > ${dataset}_${outname}_contiguator.out
	done
done

