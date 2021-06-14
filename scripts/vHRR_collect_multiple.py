import glob
import argparse
import numpy as np
import pandas as pd
from vHRR_functions import *
from network_inference import *


#command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(dest="go_file", help="input GO annotations")
parser.add_argument(dest="tmp", help="tmp folder")
parser.add_argument('-hrr', dest="hrr_sizes", help="tsv of HRR sizes per GO term")
parser.add_argument("-fdr", dest="fdr", type=float, default=0.05, help="fdr p-value cutoff")
parser.add_argument("-fc", dest="use_FC", action='store_true', help="flag to use fold change instead of p-value")            # off: False; on: True
parser.add_argument("-mqse", dest="correct_MQSE", action='store_true', help="flag to use MQSE score instead of p-value")     # off: False; on: True
args = parser.parse_args()




# PARAMETERS
sample_size=0
fdr=args.fdr
use_FC=args.use_FC
use_Fmax=True
correct_MQSE=args.correct_MQSE
use_ICu=True
extend_predictions_before_Fmax=True 
extend_predictions_after_Fmax=False
reduce_final_prediction=True


# STEP 1: READ IN STARTER PACKS
stderr.write("[INFO] Reading starter packs\n")
atlas_gene_ids = {}
all_gene_ids = []
for starter_pack in glob.glob('{}/*/starter_pack.npz'.format(args.tmp)):
    atlas = starter_pack.split('/')[-2]
    container = np.load(starter_pack)
    atlas_gene_ids[atlas] = list(container['gene_ids'])
    all_gene_ids += list(container['gene_ids'])
    del container
all_gene_ids = sorted(set(all_gene_ids))                    # be carefull with this...


# Read in GO data & calculate ICu of each GO term. Extend to be sure.
stderr.write("[INFO] Reading feature data\n")
feature_matrix_com, ICu_matrix_com, feature_gene_ids_com, go_ids_com = load_feature_matrix(args.go_file)
if extend_predictions_before_Fmax:
    extend_with_go(feature_matrix_com, go_ids_com, smaller_is_better=False)


# Match identifiers between compendia and feature matrix
stderr.write("[INFO] Matching identifiers\n")
common_ids = sorted(set(all_gene_ids) & set(feature_gene_ids_com))

if not common_ids:
    raise AssertionError('no overlap between gene ids in atlas and feature file')

# mapping dicts
compendium_ind = dict(zip(all_gene_ids, list(range(len(all_gene_ids)))))
feature_gene_ind_com = dict(zip(feature_gene_ids_com, list(range(len(feature_gene_ids_com)))))

# Remove genes from feature matrix if they are not in the compendia
feature_matrix_com = feature_matrix_com[[feature_gene_ind_com[x] for x in common_ids]]

# Add empty rows for featureless genes. Necessary to match network column indices
full_feature_matrix = np.zeros([len(all_gene_ids), len(go_ids_com)])
full_feature_matrix[[compendium_ind[x] for x in common_ids]] = feature_matrix_com



# STEP 2: GET HRR CUT-OFF VALUES
stderr.write("[INFO] Reading HRR cutoffs\n")
HRR_sizes = pd.read_csv(args.hrr_sizes, sep='\t', header=None)
HRR_sizes.columns = ['i_feature', 'GO', 'HRR', 'alpha', 'Fmax', 'atlas', 'info']
HRR_sizes = HRR_sizes.sort_values(['GO','Fmax'], ascending=True)

HRR_sizes = HRR_sizes.query("Fmax > 0")

mean_per_go = HRR_sizes.groupby('GO').mean().Fmax
HRR_sizes['select'] = HRR_sizes.apply(lambda row: True if row.Fmax >= mean_per_go[row.GO] and row.Fmax > 0 else False, axis=1)
HRR_sizes_filtered = HRR_sizes.query("select == True").set_index("GO")
mean_HRR_size = round(HRR_sizes_filtered.query("Fmax > 0")['HRR'].mean())

HRR_sizes_filtered.to_csv('{}/GO_assignment.tsv'.format(args.tmp), sep='\t')


# STEP 3: collect results

stderr.write("[INFO] Collecting individual enrichment results\n")
p_values = np.ones([len(all_gene_ids), len(go_ids_com)], dtype=np.float32)
enrichment_folds = np.ones([len(all_gene_ids), len(go_ids_com)], dtype=np.float32)
module_hits = np.ones([len(all_gene_ids), len(go_ids_com)], dtype=np.float32)


for go_term, values in HRR_sizes_filtered.iterrows():

    atlas = values['atlas']
    Fmax = values['Fmax']
    i_feature = values['i_feature']

    # init full matrices
    full_pvalue_matrix = np.ones([len(all_gene_ids), 1], dtype=np.float32)
    full_enrichments_matrix = np.ones([len(all_gene_ids), 1], dtype=np.float32)
    full_modulehits_matrix = np.ones([len(all_gene_ids), 1], dtype=np.float32)

    if Fmax > 0:
    
        enrichment_file = '{}/{}/enrichments/{}_enrichments.npz'.format(args.tmp, atlas, go_term)

        container = np.load(enrichment_file)
        i_feature = container['i_feature']
        gene_ids = atlas_gene_ids[atlas]
        assert go_ids_com[i_feature] == container['feature']
        assert i_feature == values['i_feature']
        assert go_term == container['feature']

        # fill with predicted genes
        full_pvalue_matrix[[compendium_ind[x] for x in gene_ids]] = container['p_values']
        full_enrichments_matrix[[compendium_ind[x] for x in gene_ids]] = container['enrichment_folds']
        full_modulehits_matrix[[compendium_ind[x] for x in gene_ids]] = container['module_hits']

    p_values[:, i_feature] = np.minimum(p_values[:, i_feature], full_pvalue_matrix.flatten())
    enrichment_folds[:, i_feature] = np.maximum(enrichment_folds[:, i_feature], full_enrichments_matrix.flatten())
    module_hits[:, i_feature] = np.maximum(module_hits[:, i_feature], full_modulehits_matrix.flatten())


# Run the final prediction
stderr.write("[INFO] Evaluating predictions\n")
evaluate_enrichments_results = evaluate_enrichments(p_values, 
                                                    enrichment_folds, 
                                                    full_feature_matrix, 
                                                    ICu_matrix_com, go_ids_com, 
                                                    fdr=fdr, 
                                                    use_FC=use_FC,
                                                    use_Fmax=use_Fmax, 
                                                    correct_MQSE=correct_MQSE, 
                                                    use_ICu=use_ICu, 
                                                    extend_after_predicting=extend_predictions_after_Fmax, 
                                                    reduce_final_prediction=reduce_final_prediction)

predictions, recall, precision, Fmax, max_p, predicted_per_gene, max_ic_per_gene, score_cutoff = evaluate_enrichments_results
assert predictions.shape == (len(all_gene_ids), len(go_ids_com))


avg_predicted_per_gene = np.mean(predicted_per_gene)
avg_ic_per_gene = np.mean(max_ic_per_gene[np.greater(max_ic_per_gene, 0)])


# Write predictions
stderr.write("[INFO] Writing output\n")
print("# ---- #")
print("# HRR\trecall\tprecision\tFmax\tscore_cutoff\tmax_p\tavg_predicted_per_gene\tavg_ic_per_gene")
print("# "+"\t".join(map(str, [mean_HRR_size, recall, precision, Fmax, score_cutoff, max_p, avg_predicted_per_gene, avg_ic_per_gene])))
print("# ---- #")
print("GO\tgene\tp-val\tenrichment\tmodule_hits\tMQSEscore\tICu_score\trecovered")
for i_gene in range(len(predictions)):
    for i_go in range(predictions.shape[1]):
        if not predictions[i_gene][i_go]:
            continue
        gene_id = all_gene_ids[i_gene]
        go_id = go_ids_com[i_go]
        print(go_id+"\t"+gene_id+"\t"+str(p_values[i_gene][i_go])+"\t"+str(enrichment_folds[i_gene][i_go])+"\t"+str(int(module_hits[i_gene][i_go]))+"\t"+str(np.power(10, enrichment_folds[i_gene][i_go]*np.log10(p_values[i_gene][i_go])))+"\t"+str(ICu_matrix_com[0][i_go])+"\t"+str(int(full_feature_matrix[i_gene][i_go])>0))
