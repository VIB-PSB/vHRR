#!/usr/bin/python
import go_manipulations
from os import path
from sys import argv

gene_go_file = argv[1]
filter_category = argv[2]
go_obo = argv[3]

go_tree = go_manipulations.GOtree(go_obo)

go_genes = go_manipulations.GOgenes(gene_go_file, go_tree)
go_tree.filter_category(go_genes, filter_category)
go_genes.write()
