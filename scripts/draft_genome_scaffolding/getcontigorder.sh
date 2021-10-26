#!/bin/bash
while read line; do
	strain=`echo -e "$line"| cut -f 1`
	bestmatch=`echo -e "$line"| cut -f 2`
	#../../../../${strain}_${match}Map_1/
	grep 'systematic_id' ./${strain}_${bestmatch}Map_1/PseudoContig.embl -B1 | sed 's/FT.*Contig\s\+//g' | sed 's/FT.*_id=/\t/g' | tr -d '\n' | sed 's/--/\n/g;s/"//g' | sed 's/^[0-9]\+.*\s\+/forward\t/g;s/^complement.*\s\+/reverse\t/g' > ${strain}_${bestmatch}_contigorder.out
done < bestlist 
