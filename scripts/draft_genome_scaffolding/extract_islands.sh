#!/bin/bash
while read line; do
	strain=`echo -e "$line" | cut -f 1`
	islandtype=`echo -e "$line" | cut -f 2`
	sistart=`echo -e "$line" | cut -f 3`
	siend=`echo -e "$line" | cut -f 4`
	echo "$islandtype"
	python Genbank_slicer.py -g ./${strain}_reordermerge.gbk -o ./${strain}.${islandtype}.gbk -s $sistart -e $siend 
done < SI_att_position_table 
