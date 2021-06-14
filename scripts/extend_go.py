#!/usr/bin/python
import go_manipulations
from os import path
from sys import argv

gene_go_file = argv[1]
go_obo = argv[2]

go_tree = go_manipulations.GOtree(go_obo)

gene_go = go_manipulations.GOgenes(gene_go_file, go_tree)
go_tree.extend(gene_go)
gene_go.write()

