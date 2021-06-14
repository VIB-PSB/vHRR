import argparse
import numpy as np
import pandas as pd
from vHRR_functions import *
from network_inference import *


#command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(dest="atlas", help="location of atlas or network (input)")
parser.add_argument(dest="GO_file", help="location of GO-gene-ICu file (input)")
parser.add_argument(dest="starter_pack", help="GO-wise starter pack (output)")
parser.add_argument("-gs", dest="GS", default=None, help="gold standard GO-gene file")
parser.add_argument("-sample_weighing", dest="build_network", action='store_false', help="deactivate network building, when using sample weighing")             # off: True; on: False
args = parser.parse_args()



# PARAMETERS
gold_standard_go_file = args.GS
go_file = args.GO_file


# GENERATE HRR NETWORK

# Read in the compendium
if args.atlas[-4:] == '.tsv':
    stderr.write("[INFO] Reading in compendium file\n")
    compendium, gene_ids, sample_ids = read_compendium(args.atlas, expressed_only=True)
elif args.atlas[-4:] == '.npz':
    stderr.write("[INFO] Reading HRR network\n")
    hrr_network, gene_ids = read_network_npz(args.atlas)
else:
    raise ValueError("extention of 'atlas' argument should be '.tsv' (expression atlas) or '.npz' (precomputed network))")

if args.build_network:

    if args.atlas[-4:] == '.tsv':

        # Turn compendium into PCC networks
        stderr.write("[INFO] Calculating PCC networks\n")
        pcc_network = simple_correlation(compendium)
        del compendium

        # Convert to HRR networks
        stderr.write("[INFO] Converting to HRR networks\n")
        hrr_network = convert_hrr(pcc_network)
        del pcc_network

    full_hrr_network = np.copy(hrr_network)


# Read in GO data & calculate ICu of each GO term. Extend to be sure.
stderr.write("[INFO] Reading feature data\n")
feature_matrix_com, ICu_matrix_com, feature_gene_ids_com, go_ids_com = load_feature_matrix(go_file)
extend_with_go(feature_matrix_com, go_ids_com, smaller_is_better=False)


# Match identifiers between compendia and feature matrix
stderr.write("[INFO] Matching identifiers\n")
common_ids = sorted(set(gene_ids) & set(feature_gene_ids_com))

if not common_ids:
    raise AssertionError('no overlap between gene ids in atlas and feature file')

# mapping dicts
compendium_ind = dict(zip(gene_ids, list(range(len(gene_ids)))))
feature_gene_ind_com = dict(zip(feature_gene_ids_com, list(range(len(feature_gene_ids_com)))))

# select relevant network rows (ie genes that have a input GO annotation; all cols are retained, in order to evaluate using the true clusters)
if args.build_network:
    hrr_network = hrr_network[[compendium_ind[x] for x in common_ids]]
    if not hrr_network.any():
        raise ValueError('HRR network broken')
else:
    hrr_network = None
    full_hrr_network = None

# Remove genes from feature matrix if they are not in the compendium
feature_matrix_com = feature_matrix_com[[feature_gene_ind_com[x] for x in common_ids]]

# Add empty rows for featureless genes. Necessary to match network column indices
full_feature_matrix = np.zeros([len(gene_ids), len(go_ids_com)])
full_feature_matrix[[compendium_ind[x] for x in common_ids]] = feature_matrix_com

stderr.write("[INFO] Preparing gold standard\n")
if gold_standard_go_file == None:
    selected_gene_ind = list(range(len(common_ids)))
    selected_gene_ids = np.array(common_ids)
    gold_standard = feature_matrix_com[selected_gene_ind]
else:
    selected_gene_ids = np.array(common_ids)
    gold_standard = read_gold_standard_from_file(gold_standard_go_file, selected_gene_ids, go_ids_com)


stderr.write("[INFO] Saving starter pack\n")
np.savez(args.starter_pack, 
         go_ids_com=go_ids_com,
         hrr_network_copy_start=hrr_network,
         full_hrr_network = full_hrr_network,
         full_feature_matrix=full_feature_matrix, 
         selected_gene_ids=selected_gene_ids, 
         gene_ids=gene_ids, 
         common_ids=np.array(common_ids),
         gold_standard=gold_standard,
         ICu_matrix_com=ICu_matrix_com)

for el in go_ids_com:
    print(el)