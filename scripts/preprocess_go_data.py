import os
import pandas as pd 
from sys import argv, stderr
from backend import get_obo


def extend_annotations(go_file, out_file, go_obo):

    # extend predictions
    os.system("scripts/extend_go.py {} {} > tmp/intermediate_go_file1.tmp".format(go_file, go_obo))

    # filter BP terms
    os.system("scripts/filter_go_category.py tmp/intermediate_go_file1.tmp biological_process {} > tmp/intermediate_go_file2.tmp".format(go_obo))

    # drop 'biological_process'
    os.system("grep -v 'GO:0008150' tmp/intermediate_go_file2.tmp > tmp/intermediate_go_file3.tmp")

    # cut 3rd column
    os.system("cut -f -2 tmp/intermediate_go_file3.tmp > {}".format(out_file))

    # delete tmp files
    os.system("rm tmp/intermediate_go_file*.tmp")


def calculate_IC(go_file, out_file):
    os.system("scripts/calculate_ICu_per_GO.py {} > {}".format(go_file, out_file))


def addIC(go_gene, go_ic, out_file):
    
    # A
    go_gene = pd.read_csv(go_gene, sep='\t', header=None)
    go_gene.columns = ['GO', 'gene_id']
    go_gene = go_gene.set_index('GO')
    
    # B
    go_ic = pd.read_csv(go_ic, sep='\t', header=None)
    go_ic.columns = ['GO', 'icu']
    go_ic = go_ic.set_index('GO')
    
    # join B to A
    go_gene_ic = go_gene.join(go_ic, how='left')
    assert len(go_gene) == len(go_gene_ic)
    
    # write
    go_gene_ic.to_csv(out_file, sep='\t', header=None)




if __name__ == '__main__':


    go_file = argv[1]
    go_gene_ic = argv[2]

    go_obo = get_obo()

    # make temp dir
    os.system("mkdir -p tmp")

    # extend predictions + filter BP
    stderr.write('[INFO]\tExtending GO-gene\n')
    extend_annotations(go_file, 'tmp/extended_go_file.tmp', go_obo)

    # calculate IC
    stderr.write('[INFO]\tCalculating GO information content\n')
    calculate_IC('tmp/extended_go_file.tmp', 'tmp/ic_per_go.tsv')

    # join dfs
    stderr.write('[INFO]\tAdding GO information content column\n')
    addIC('tmp/extended_go_file.tmp', 'tmp/ic_per_go.tsv', go_gene_ic)

    os.system("rm tmp/extended_go_file.tmp")

