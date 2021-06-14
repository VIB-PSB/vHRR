import glob
import argparse
import numpy as np
import pandas as pd
from vHRR_functions import *
from network_inference import *


#command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(dest="atlas", help="location of atlas or network")
parser.add_argument(dest="starter_pack", help="GO-wise starter pack")
parser.add_argument('-hrr', dest="hrr_sizes", help="tsv of HRR sizes per GO term")
parser.add_argument('-g', dest="global_hrr", default=None, help="global optimal HRR. By default, the mean of GO-specific HRR values is taken")
parser.add_argument("-fdr", dest="fdr", type=float, default=0.05, help="fdr p-value cutoff")
parser.add_argument("-fc", dest="use_FC", action='store_true', help="flag to use fold change instead of MQSE score")            # off: False; on: True
parser.add_argument("-mqse", dest="correct_MQSE", action='store_true', help="flag to use MQSE score")                           # off: False; on: True
parser.add_argument('-enrichments', dest="enrichments", default=None, help="folder with enrichment npzs per GO term")
parser.add_argument("-all_features", dest="all_features", action='store_true', help="do predictions for all GO terms (not only terms where TPs were found)")             # off: False; on: True
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


# STEP 1: READ IN STARTER PACK
container = np.load(args.starter_pack)
go_ids_com = container['go_ids_com']
full_feature_matrix = container['full_feature_matrix']
ICu_matrix_com = container['ICu_matrix_com']
common_ids = list(container['common_ids'])
del container


# Read in the compendium
#stderr.write("[INFO] Reading in compendium file\n")
#compendium, gene_ids, sample_ids = read_compendium(args.atlas, expressed_only=True)
if args.atlas[-4:] == '.tsv':
    stderr.write("[INFO] Reading in compendium file\n")
    compendium, gene_ids, sample_ids = read_compendium(args.atlas, expressed_only=True)
elif args.atlas[-4:] == '.npz':
    stderr.write("[INFO] Reading HRR network\n")
    full_hrr_network, gene_ids = read_network_npz(args.atlas)
else:
    assert False
compendium_ind = dict(zip(gene_ids, list(range(len(gene_ids)))))


# STEP 2: GET HRR CUT-OFF VALUES
if args.hrr_sizes:
    HRR_sizes = pd.read_csv(args.hrr_sizes, sep='\t', header=None)
    HRR_sizes.columns = ['i_feature', 'feature', 'HRR', 'alpha', 'Fmax', 'atlas', 'info']
    HRR_sizes = HRR_sizes.sort_values(by='i_feature')
    mean_HRR_size = round(HRR_sizes.query("HRR > 0")['HRR'].mean())
    stderr.write('Mean HRR cut-off: {}\n'.format(np.mean(mean_HRR_size)))

    if args.global_hrr:
        global_hrr = int(args.global_hrr)
    else:
        global_hrr = mean_HRR_size
    recovered_terms = set(HRR_sizes.query("HRR > 0").i_feature)
    HRR_sizes.loc[HRR_sizes.HRR == 0, 'HRR'] = global_hrr
    HRR_sizes = HRR_sizes['HRR'].values
else:
    assert args.enrichments
    mean_HRR_size = np.nan


if args.enrichments:

    p_values = np.ones([len(gene_ids), len(go_ids_com)], dtype=np.float32)
    enrichment_folds = np.ones([len(gene_ids), len(go_ids_com)], dtype=np.float32)
    module_hits = np.ones([len(gene_ids), len(go_ids_com)], dtype=np.float32)

    for enrichment_file in glob.glob('{}/*_enrichments.npz'.format(args.enrichments)):

        container = np.load(enrichment_file)
        i_feature = container['i_feature']
        assert go_ids_com[i_feature] == container['feature']

        p_values[:, i_feature] = container['p_values'].flatten()
        enrichment_folds[:, i_feature] = container['enrichment_folds'].flatten()
        module_hits[:, i_feature] = container['module_hits'].flatten()

    # add unannotated genes
    #p_values = np.ones([len(gene_ids), len(go_ids_com)])
    #p_values[[compendium_ind[x] for x in common_ids]] = p_values_small

    #enrichment_folds = np.ones([len(gene_ids), len(go_ids_com)])
    #enrichment_folds[[compendium_ind[x] for x in common_ids]] = enrichment_folds_small



if not args.enrichments or args.all_features:
    # GENERATE HRR NETWORK


    if args.atlas[-4:] == '.tsv':
        # Turn compendium into PCC networks
        stderr.write("[INFO] Calculating PCC networks\n")
        pcc_network = simple_correlation(compendium)
        del compendium

        # Convert to HRR networks
        stderr.write("[INFO] Converting to HRR networks\n")
        full_hrr_network = convert_hrr(pcc_network)
        del pcc_network


    # STEP 3: CALCULATE PREDICTIONS

    # cutoff network at mean or global HRR
    global_network = np.copy(full_hrr_network)
    cutoff_network(global_network, global_hrr)

    # Get p-values and enrichment folds
    stderr.write("[INFO] Calculating enrichments\n")
    new_p_values, new_enrichment_folds, new_module_hits = calculate_enrichments_final(full_hrr_network, 
                                                                                      global_network,
                                                                                      full_feature_matrix, 
                                                                                      HRR_sizes, 
                                                                                      global_hrr,
                                                                                      recovered_terms,
                                                                                      gene_ids, 
                                                                                      gene_ids,
                                                                                      only_global=args.all_features)

    if args.all_features:

        #test_array = p_values[np.equal(p_values, 1)]+new_p_values[np.equal(new_p_values, 1)]
        #test_array = p_values+new_p_values
        #assert np.all(test_array[np.greater_equal(test_array, 1)])

        p_values *= new_p_values
        enrichment_folds *= new_enrichment_folds
        module_hits *= new_module_hits
    else:
        p_values = new_p_values
        enrichment_folds = new_enrichment_folds
        module_hits = new_module_hits


    # Extend p-values to parental terms if required
    stderr.write("[INFO] Extending enrichments\n")
    if extend_predictions_before_Fmax:
        extend_with_go(p_values, go_ids_com, smaller_is_better=True)
        extend_with_go(enrichment_folds, go_ids_com, smaller_is_better=False)
        extend_with_go(module_hits, go_ids_com, smaller_is_better=False)


    #np.savez('arabidopsis/sample_weighing/gene1-500_enrichments_parallel.npz', 
    #         p_values=p_values,
    #         enrichment_folds=enrichment_folds)
    #assert False

# Run the final prediction
stderr.write("[INFO] Generating predictions\n")
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
        gene_id = gene_ids[i_gene]
        go_id = go_ids_com[i_go]
        print(go_id+"\t"+gene_id+"\t"+str(p_values[i_gene][i_go])+"\t"+str(enrichment_folds[i_gene][i_go])+"\t"+str(int(module_hits[i_gene][i_go]))+"\t"+str(np.power(10, enrichment_folds[i_gene][i_go]*np.log10(p_values[i_gene][i_go])))+"\t"+str(ICu_matrix_com[0][i_go])+"\t"+str(int(full_feature_matrix[i_gene][i_go])>0))
