#!/usr/bin/env python
# Script by JK: https://github.com/kwongj/sort-contigs/blob/master/sort-contigs.py
#modified to also look at direction
# Sort contigs in GBK or FASTA file

# Use modern print function from python 3.x
from __future__ import print_function

# Import modules
import argparse
import os
import sys
import StringIO
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq, Alphabet
# Usage
parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Sort contigs in a GBK or FASTA file',
		usage='\n  %(prog)s [--order order.txt] CONTIGS.GBK | CONTIGS.FA')
parser.add_argument('contigs', metavar='FILE', nargs=1, help='genbank or FASTA file to sort')
parser.add_argument('--fmt', metavar='FORMAT', nargs=1, required=True, help='file format (fasta|genbank)')
parser.add_argument('--order', metavar='FILE', nargs=1, help='specified order (list of contigs must match contig names)')
parser.add_argument('--out', metavar='FILE', nargs=1, help='Output file (optional - otherwise will print to stdout)')
parser.add_argument('--version', action='version', version='%(prog)s version 1.0\nUpdated 23-Nov-2016\nScript by JK')
args = parser.parse_args()

ctgfile = args.contigs[0]
format = args.fmt[0]
if args.order:
	orderfile = args.order[0]
if args.out:
	outfile = args.out[0]

# Functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Check file exists
def check_file(f):
	if os.path.isfile(f) == True:
		err('ERROR: Output file "{}" already exists. Please specify a new output file.'.format(f))

# Parse contigs
ctglist = {}
ctg_ids = []
seqNEW = []
dirdict = {}
for seqREC in SeqIO.parse(ctgfile, format):
	ctglist[seqREC.id] = seqREC
	ctg_ids.append(seqREC.id)
if len(ctg_ids) < 1:
	err('ERROR: Check file format.')
if args.order:
	ctg_ids=[]
	with open(orderfile) as file:
		for line in file.read().splitlines():
			tmp=line.split('\t')
			ctg_ids.append(tmp[1])
			dirdict[tmp[1]] = tmp[0]
			#ctg_ids = [line.rstrip() for line in file]
else:
	ctg_ids.sort()
for id in ctg_ids:
	if id in dirdict.keys():
		if dirdict[id] == "reverse":
			seqNEW.append(ctglist[id].reverse_complement(id=True,description=True,annotations=True,dbxrefs=True))
		else:
			seqNEW.append(ctglist[id])

#merge genbank seqNEW into single record
new_rec = seqNEW[0]

spacer_seq = UnknownSeq(1000, alphabet=Alphabet.generic_dna, character='N')
rec_id = new_rec.id
rec_name = new_rec.name
rec_desc = new_rec.description
date = new_rec.annotations.get('date', '')
source = new_rec.annotations.get("source", '')
organism = new_rec.annotations.get('organism', '')
taxonomy = new_rec.annotations.get('taxonomy', [])
data_file_division = new_rec.annotations.get('data_file_division', 'UNK')
topology = new_rec.annotations.get('topology', 'linear')

for i, rec in enumerate(seqNEW[1:]):
	spacer_id = 'spacer_{}'.format(i + 1)

	spacer_feature = SeqFeature(FeatureLocation(0, 1000, 0),
								type='misc_feature', id=spacer_id,
								qualifiers={'note': [spacer_id]})

	spacer_rec = SeqRecord(spacer_seq, id=spacer_id, name=spacer_id,
						   description=spacer_id, features=[spacer_feature])

	new_rec = new_rec + spacer_rec + rec

new_rec.id = rec_id
new_rec.name = rec_name
new_rec.description = rec_desc
new_rec.annotations["date"] = date
new_rec.annotations["source"] = source
new_rec.annotations["organism"] = organism
new_rec.annotations["taxonomy"] = taxonomy
new_rec.annotations["data_file_division"] = data_file_division
new_rec.annotations["topology"] = topology

# Write new ordered file or print to stdout
if args.out:
	check_file(outfile)
	msg('Sorted contigs saved to "{}" ... '.format(outfile))
	SeqIO.write(new_rec, outfile, format)
	#SeqIO.write(seqNEW, outfile, format)
else:
	seqFILE = StringIO.StringIO()
	SeqIO.write(new_rec, seqFILE, format)
	#SeqIO.write(seqNEW, seqFILE, format)
	output = seqFILE.getvalue().rstrip()
	print(output)

sys.exit(0)
