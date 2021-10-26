#!/bin/bash
for file in `ls -1 *contigorder.out`; do
	strain=`echo -e "$file" | sed 's/Map_chromosome//g;s/_.*//g'`
	#strain=`head -n1 $file | cut -f 2 | sed 's/_[0-9]\+$//g'`
	python combinecontigsorder.py --fmt genbank --order $file ./draft_db/${strain}.gbk > ${strain}_reordermerge.gbk
done
