module load python

INFILE=input/go_gene.tsv
OUTFILE=input/go_gene_ic.tsv

python3 scripts/preprocess_go_data.py $INFILE $OUTFILE
