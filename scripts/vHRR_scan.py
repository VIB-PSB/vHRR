import argparse
import numpy as np
import pandas as pd
import stamakro_mlc as mlc
from vHRR_functions import *
from network_inference import *

#command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(dest="atlas", help="location of atlas or network (input)")
parser.add_argument(dest="starter_pack", help="GO-wise starter pack")
parser.add_argument(dest="go_term", help="current GO term")
parser.add_argument("-r", dest="hrr_range", default='10,401,50', help="hrr range to be scanned: <min>,<max>,<initial_step>")
parser.add_argument("-hrr", dest="HRR", type=int, default=None, help="predefined HRR threshold")
parser.add_argument("-fdr", dest="fdr", type=float, default=0.05, help="fdr p-value cutoff")
parser.add_argument("-fc", dest="use_FC", action='store_true', help="flag to use fold change instead of MQSE score") 
parser.add_argument("-mqse", dest="correct_MQSE", action='store_true', help="flag to use p-values instead of MQSE score")       
parser.add_argument("-sample_weighing", dest="sample_weighing", action='store_true', help="activate sample weighing")
parser.add_argument("-out_folder", dest="out_folder", default=None, help="output folder for enrichment npz")
args = parser.parse_args()



# PARAMETERS
hrr_range=list(range(*[int(el) for el in args.hrr_range.split(',')]))
sample_size=0
fdr=args.fdr
use_FC=args.use_FC
use_Fmax=True
correct_MQSE=args.correct_MQSE
use_ICu=False
extend_predictions_before_Fmax=True 
extend_predictions_after_Fmax=False
reduce_final_prediction=True
randomize_network = False
randomize_GO = False



# STEP 1: READ IN STARTER PACK
container = np.load(args.starter_pack)
go_ids_com = container['go_ids_com']
hrr_network_copy_start = container['hrr_network_copy_start']
full_feature_matrix = container['full_feature_matrix']
selected_gene_ids = container['selected_gene_ids'] 
common_ids = list(container['common_ids'])
gold_standard = container['gold_standard']
ICu_matrix_com = container['ICu_matrix_com']
del container



# STEP 2: READ IN GO TERM

feature_indexes = dict(zip(go_ids_com, list(range(len(go_ids_com)))))

feature = args.go_term
i_feature = feature_indexes[feature]

stderr.write('Feature {} {}\n'.format(i_feature, feature))

# number of genes for this GO
feature_size_here = np.sum(full_feature_matrix[:,i_feature])

if feature_size_here < 2:
    stderr.write('[INFO] feature size = 1 -> skipping term {}\n'.format(feature_size_here))    
    print('{}\t{}\t{}\t{}\t{}\t{}\t[INFO] feature size = 1 -> skipping term'.format(i_feature, feature, 0, 0, 0, args.atlas.split('/')[-1].split('.')[0]))
    exit()

stderr.write('Number of feature genes: {}\n'.format(feature_size_here))



alpha_scores = {}
alpha_size = {}
weights_dict = {}

# Loop hrr sizes from large to small    
for alpha in [i/10 for i in range(0,11)]:

    stderr.write('\n[INFO] Setting alpha at {}\n'.format(alpha))
    hrr_range = list(range(*[int(el) for el in args.hrr_range.split(',')]))
    skip_alpha = False

    # GENERATE HRR NETWORK
    #if args.build_network:
    if args.sample_weighing:

        try:
            del hrr_network
            del hrr_network_copy_start
        except:
            pass

        # Read in the compendium
        stderr.write("[INFO] Reading in compendium file\n")
        compendium, gene_ids, sample_ids = read_compendium(args.atlas, expressed_only=True)
        compendium_ind = dict(zip(gene_ids, list(range(len(gene_ids)))))
        

        #if args.sample_weighing:
        
        # sample weighing    
        stderr.write("[INFO] Weighing samples\n")
        sample_weights = mlc.mlcTrain(compendium, full_feature_matrix[:,i_feature], a=alpha)
        if sample_weights.max() == 0:
            stderr.write("[INFO] alpha {} skipped, all weights 0\n\n".format(alpha))
            continue
        sample_weights = sample_weights/sample_weights.max()
        sample_weights[np.equal(sample_weights, 0)] = 1/1000000
        
        for other_alpha, other_weights in weights_dict.items():
            weights_pcc = simple_correlation(np.array([sample_weights, other_weights]))
            
            if weights_pcc[1][0] > 0.99:
                stderr.write("[INFO] alpha {} skipped, equal to alpha {} (PCC = {})\n\n".format(alpha, other_alpha, weights_pcc[1][0]))
                skip_alpha = True
                break

        weights_dict[alpha] = sample_weights

        if skip_alpha:
            continue

        # Turn compendium into PCC networks
        stderr.write("[INFO] Calculating PCC networks\n")
        pcc_network = weighted_correlation(compendium, sample_weights).astype(np.float32)
        #else:
        #    stderr.write("[INFO] Calculating PCC networks\n")
        #    pcc_network = simple_correlation(compendium)
        del compendium

        
        # Convert to HRR networks
        stderr.write("[INFO] Converting to HRR networks\n")
        hrr_network = convert_hrr(pcc_network)
        del pcc_network

        # Randomize this network if requested
        if randomize_network:
            for i in range(len(hrr_network)):
                np.random.shuffle(hrr_network[i])


        # select relevant network rows (ie genes that have a input GO annotation; all cols are retained, in order to evaluate using the true clusters)
        hrr_network = hrr_network[[compendium_ind[x] for x in common_ids]]
        hrr_network_copy_start = np.copy(hrr_network)
            

    else:
        hrr_network = np.copy(hrr_network_copy_start)

    if not hrr_network.any() or np.isnan(hrr_network).all():
        stderr.write('\n[INFO] HRR network broken, skipping alpha {}\n'.format(alpha))
        continue


    # STEP 3: SCAN FIRST RANGE OF HRR SIZES

    hrr_scores = dict() # { HRR_size : Fmax }

    # Loop hrr sizes from large to small    
    for hrr_size in sorted(hrr_range, reverse=True):     

        # Set the hrr_size cutoff on the hrr networks
        stderr.write('\n[INFO] Setting HRR size at {}\n'.format(hrr_size))
        cutoff_network(hrr_network, hrr_size)

        # Calculate enrichments for true compendium
        stderr.write('[INFO] Generating enrichments\n')
        p_values, enrichment_folds, module_hits = calculate_enrichments_scan(hrr_network, full_feature_matrix[:,i_feature], i_feature, common_ids, selected_gene_ids)

        if p_values.min() > 0.05:
            hrr_scores[hrr_size] = 0
            stderr.write('[INFO] all p-values = 1, skipping HRR cut-off\n')
            continue

        # Extend p-values to parental terms if required
        if extend_predictions_before_Fmax:
            extend_with_go(p_values, [feature], smaller_is_better=True)
            extend_with_go(enrichment_folds, [feature], smaller_is_better=False)
            extend_with_go(module_hits, [feature], smaller_is_better=False)


        # Evaluate predictions with the FDR cutoff
        stderr.write('[INFO] Evaluating predictions\n')
        evaluate_enrichments_results = evaluate_enrichments(p_values, 
                                                            enrichment_folds, 
                                                            gold_standard[:, np.newaxis, i_feature], 
                                                            ICu_matrix_com[np.newaxis, :, i_feature], 
                                                            [feature], 
                                                            scanning=True, fdr=fdr, use_FC=use_FC, use_Fmax=use_Fmax, 
                                                            correct_MQSE=correct_MQSE, use_ICu=use_ICu, 
                                                            extend_after_predicting=extend_predictions_after_Fmax, 
                                                            reduce_final_prediction=reduce_final_prediction)

        predictions, recall, precision, Fmax, max_p, predicted_per_gene, max_ic_per_gene, score_cutoff = evaluate_enrichments_results


        hrr_scores[hrr_size] = Fmax
        stderr.write('[INFO] Fmax for HRR{} = {} (R: {} -- P: {}) > cutoff = {}\n'.format(hrr_size, Fmax, recall, precision, score_cutoff))


    if max(hrr_scores.values()) == 0:
        if args.sample_weighing:
            stderr.write('[INFO] No TPs, skipping alpha\n')
            continue
        else:
            stderr.write('[INFO] No TPs, skipping term\n')
            print('{}\t{}\t{}\t{}\t{}\t{}\t[INFO] No TPs, skipping term'.format(i_feature, feature, 0, 0, 0, args.atlas.split('/')[-1].split('.')[0]))
            exit()
    elif max(hrr_scores.values()) == min(hrr_scores.values()) and args.sample_weighing:
        stderr.write('[INFO] No difference in Fmax, skipping alpha\n')
        continue



    # STEP 4: ZOOM IN ON OPTIMAL HRR

    if len(hrr_range) > 1:
        range_start = min(hrr_range)
        range_stop = max(hrr_range)
        step = hrr_range[1] - hrr_range[0]
        stderr.write('[INFO] Zooming in on optimal HRR\n')

        while step > 1:
            step = max(1, int(step/2))
            best_so_far = sorted(hrr_scores.keys(), key=lambda x: hrr_scores[x])[-1]
            next_steps = sorted([x for x in [best_so_far+step, best_so_far-step] if range_start < x and x < range_stop], reverse=True)
            hrr_network = np.copy(hrr_network_copy_start)

            for hrr_size in next_steps:

                # Cut off network here
                stderr.write('\n[INFO] Setting HRR size at {}\n'.format(hrr_size))
                cutoff_network(hrr_network, hrr_size)

                # Calculate enrichments for true compendium
                stderr.write('[INFO] Generating enrichments\n')
                p_values, enrichment_folds, module_hits = calculate_enrichments_scan(hrr_network, full_feature_matrix[:,i_feature], i_feature, common_ids, selected_gene_ids)

                # Extend p-values to parental terms if required
                if extend_predictions_before_Fmax:
                    extend_with_go(p_values, [feature], smaller_is_better=True)
                    extend_with_go(enrichment_folds, [feature], smaller_is_better=False)
                    extend_with_go(module_hits, [feature], smaller_is_better=False)

                # Evaluate predictions with the FDR cutoff
                stderr.write('[INFO] Evaluating predictions\n')
                evaluate_enrichments_results = evaluate_enrichments(p_values, 
                                                                    enrichment_folds, 
                                                                    gold_standard[:, np.newaxis, i_feature], 
                                                                    ICu_matrix_com[np.newaxis, :, i_feature], 
                                                                    [feature], 
                                                                    scanning=True, fdr=fdr, use_FC=use_FC, use_Fmax=use_Fmax, 
                                                                    correct_MQSE=correct_MQSE, use_ICu=use_ICu, 
                                                                    extend_after_predicting=extend_predictions_after_Fmax, 
                                                                    reduce_final_prediction=reduce_final_prediction)

                predictions, recall, precision, Fmax, max_p, predicted_per_gene, max_ic_per_gene, score_cutoff = evaluate_enrichments_results


                hrr_scores[hrr_size] = Fmax

                stderr.write('[INFO] Fmax for HRR{} = {} (R: {} -- P: {}) > cutoff = {}\n'.format(hrr_size, Fmax, recall, precision, score_cutoff))
                #stderr.flush()

    hrr_size = sorted(hrr_scores.keys(), key=lambda x: hrr_scores[x])[-1]
    stderr.write('[INFO] Found optimal HRR at {} for alpha {}\n\n'.format(hrr_size, alpha))
    
    alpha_scores[alpha] = hrr_scores[hrr_size]
    alpha_size[alpha] = hrr_size

    if not args.sample_weighing:
        break


if not alpha_scores:
    stderr.write('[INFO] No TPs, skipping term\n')
    print('{}\t{}\t{}\t{}\t{}\t{}\t[INFO] No TPs, skipping term'.format(i_feature, feature, 0, 0, 0, args.atlas.split('/')[-1].split('.')[0]))
    exit()

alpha = sorted(alpha_scores.keys(), key=lambda x: alpha_scores[x])[-1]
hrr_size = alpha_size[alpha]
Fmax = alpha_scores[alpha]

assert max(alpha_scores.values()) == alpha_scores[alpha]

if max(alpha_scores.values()) == 0:
    stderr.write('[INFO] No TPs, skipping term\n')
    print('{}\t{}\t{}\t{}\t{}\t{}\t[INFO] No TPs, skipping term'.format(i_feature, feature, 0, 0, 0, args.atlas.split('/')[-1].split('.')[0]))
    exit()




stderr.write('[INFO] Found optimal alpha: {} with HRR {}\n\n'.format(alpha, hrr_size))




if args.sample_weighing:

    # Read in the compendium
    stderr.write("[INFO] Reading in compendium file\n")
    compendium, gene_ids, sample_ids = read_compendium(args.atlas, expressed_only=True)

    # sample weighing    
    stderr.write("[INFO] Weighing samples with alpha = {}\n".format(alpha))
    sample_weights = mlc.mlcTrain(compendium, full_feature_matrix[:,i_feature], a=alpha)
    sample_weights = sample_weights/sample_weights.max()
    sample_weights[np.equal(sample_weights, 0)] = 1/1000000
    
    # Turn compendium into PCC networks
    stderr.write("[INFO] Calculating PCC networks\n")
    pcc_network = weighted_correlation(compendium, sample_weights).astype(np.float32)
    del compendium

    # Convert to HRR networks
    stderr.write("[INFO] Converting to HRR networks\n")
    hrr_network = convert_hrr(pcc_network)
    del pcc_network

else:
    # Read full network
    container = np.load(args.starter_pack)
    hrr_network = container['full_hrr_network']
    gene_ids = container['gene_ids']
    del container


if not hrr_network.any():
    raise ValueError('HRR network broken')

# Cut off network
stderr.write('\n[INFO] Setting HRR size at optimal HRR: {}\n'.format(hrr_size))
cutoff_network(hrr_network, hrr_size)

# Calculate enrichments for true compendium
stderr.write('[INFO] Generating enrichments\n')
p_values, enrichment_folds, module_hits = calculate_enrichments_scan(hrr_network, full_feature_matrix[:,i_feature], i_feature, gene_ids, gene_ids)

# Extend p-values to parental terms if required
if extend_predictions_before_Fmax:
    extend_with_go(p_values, [feature], smaller_is_better=True)
    extend_with_go(enrichment_folds, [feature], smaller_is_better=False)
    extend_with_go(module_hits, [feature], smaller_is_better=False)

np.savez('{}/{}_enrichments.npz'.format(args.out_folder, feature), 
         p_values=p_values,
         enrichment_folds=enrichment_folds,
         module_hits=module_hits,
         feature=feature,
         i_feature=i_feature,
         gene_ids=gene_ids)


print('{}\t{}\t{}\t{}\t{}\t{}\t[INFO] Found optimal HRR and alpha'.format(i_feature, feature, hrr_size, alpha, Fmax, args.atlas.split('/')[-1].split('.')[0]))
