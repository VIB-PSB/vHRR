#!/usr/bin/python
from sys import argv
from math import log

go_files = argv[1:]
go_data = dict() # { go_id: set(gene_id) }
total_genes = set()
def read_go_file(go_file):
	with open(go_file, 'r') as reader:
		for line in reader:
			go_id, gene_id = line.strip().split("\t")[:2]
			if not go_id in go_data:
				go_data[go_id] = set()
			go_data[go_id].add(gene_id)
			total_genes.add(gene_id)
for go_file in go_files:
	read_go_file(go_file)

total_genes = float(len(total_genes))
for go_id in sorted(go_data):
	ICu = -1 * log(float(len(go_data[go_id])) / total_genes, 10) / log(total_genes, 10)
	print(go_id+"\t"+str(ICu))

