import numpy as np
from sys import stderr
import gc
import pdb
from math import log
from py_hyper import Hyper # Interface to dreec's C++ enricher code
from go_manipulations import GOtree
from network_inference import *
from backend import get_obo



go_obo = get_obo()



def load_feature_matrix(go_file, ICu_col=True):
    
    '''
    Reads in a gene-GO file, augmented by adding the ICu weights as an extra column.
    Outputs a numpy matrix of 1 and 0 values indicating which genes (rows) are annotated with which feature (columns).
    '''

    go_data = dict()     # { gene_id: { go_id : is_exp } }
    ICu_weights = dict() # { go_id: ICu }

    with open(go_file, 'r') as reader:
        
        # Read the file into the go_data and ICu_weights dictionaries
        for line in reader:
            
            # Skip comment lines
            if line[0] == "#":
                continue
            
            split = line.strip().split("\t")
            # Sanity check
            if ICu_col and len(split) != 3:
                stderr.write("[ERR] Bad format in GO file: Expected 3 columns but found "+str(len(split))+"\n")
                exit()
            elif not ICu_col and len(split) != 2:
                stderr.write("[ERR] Bad format in GO file: Expected 2 columns but found "+str(len(split))+"\n")
                exit()
            
            go_id = split[0]
            gene_id = split[1]
            if ICu_col: 
                ICu = float(split[2])
            else:
                ICu = 1

            # Store the data
            if not gene_id in go_data:
                go_data[gene_id] = set()
            go_data[gene_id].add(go_id)
            ICu_weights[go_id] = ICu
    
    # Filter experimental/computational terms | thats not what this block is doing
    gene_ids = sorted(go_data.keys())
    go_ids = sorted(ICu_weights.keys())
    ICu_weights = np.array([[ICu_weights[x] for x in go_ids]], dtype=np.float32)
    

    # Create the feature matrix
    feature_matrix = np.zeros((len(gene_ids), len(go_ids)), dtype=np.int32)
    for i_gene in range(len(gene_ids)):
        feature_matrix[i_gene][[x in go_data[gene_ids[i_gene]] for x in go_ids]] += 1
    
    # filter single gene terms
    #min_hits = None
    #if min_hits == 2:
    #    multiple_gene_go_indexes = np.greater(feature_matrix_com.sum(axis=0), 1)
    #    go_ids = go_ids[multiple_gene_go_indexes]
    #    ICu_matrix_com = np.array([ICu_matrix_com[0][multiple_gene_go_indexes]])
    #    feature_matrix = feature_matrix[:,multiple_gene_go_indexes]

    return feature_matrix, ICu_weights, gene_ids, go_ids




def extend_with_go(array, go_ids, smaller_is_better=True):
    
    '''
    Extends a gene x go array by copying better values from child to parent terms.
    For p-values, smaller is better. For enrichment folds, larger is better.
    '''
    
    go_tree = GOtree(go_obo)

    go_set = set(go_ids)
    go_ind = dict(zip(go_ids, list(range(len(go_ids)))))
    for i_go, go_id in enumerate(go_ids):
        ancestors = go_set & go_tree.get_ancestors(go_id)
        for ancestor in ancestors:
            i_anc = go_ind[ancestor]
            if smaller_is_better:
                array[:,i_anc] = np.minimum(array[:,i_anc], array[:,i_go])
            else:
                array[:,i_anc] = np.maximum(array[:,i_anc], array[:,i_go])




def reduce_with_go(array, go_ids):
    
    '''
    Removes parental GO terms from a prediction
    '''
    
    go_tree = GOtree(go_obo)

    go_set = set(go_ids)
    go_ind = dict(zip(go_ids, list(range(len(go_ids)))))
    for i_go, go_id in enumerate(go_ids):
        go_vector = np.greater(array[:,i_go], 0)
        ancestors = go_set & go_tree.get_ancestors(go_id)
        for ancestor in ancestors:
            i_anc = go_ind[ancestor]
            array[go_vector,i_anc] = 0


def calculate_enrichments_scan(network, feature_matrix, i_feature, gene_ids, genes_to_predict, min_hits=2):
    

    gene_ind = dict(zip(gene_ids, list(range(len(gene_ids)))))
    total_size = network.shape[1]
    

    # number of genes for this GO
    feature_size = np.sum(feature_matrix)
    

    # Prepare emtpy matrices for p-values and enrichment folds
    p_values = np.ones([len(genes_to_predict), 1], dtype=np.float32)
    enrichment_folds = np.ones([len(genes_to_predict), 1], dtype=np.float32)
    module_hits = np.ones([len(genes_to_predict), 1], dtype=np.float32)


    hyper = Hyper()

    # Loop genes
    for i_gene, gene in enumerate(genes_to_predict):
        
        # translate gene id into index (of common_ids, not the original compendium index!).
        ind = gene_ind[gene]                            
        

        # Get indexes of cluster genes
        cluster = np.logical_not(np.isnan(network[ind]))                # logical_not(isnan()): element-wise NaN -> False
        
        sample_size = np.sum(cluster)
        if sample_size == 0:
            raise ValueError('[ERROR] Empty cluster. feature {}'.format(i_feature))
        elif sample_size == total_size:
            raise ValueError('[ERROR] cluster = all genes. feature {}'.format(i_feature))

        # test if indices still match
        assert len(cluster) == len(feature_matrix)


        # number of cluster genes for this GO
        cluster_feature_hits = np.sum(feature_matrix[cluster])
        assert isinstance(cluster_feature_hits, np.float64)
        

        # skip feature when too few hits
        if cluster_feature_hits < min_hits: 
            continue
            
            
        # Get enrichment fold
        total_hits = feature_size
        sample_hits = cluster_feature_hits     
        
        
        enrichment_fold = (float(sample_hits)/sample_size) / (float(total_hits)/total_size)
        enrichment_folds[i_gene][0] = enrichment_fold
        module_hits[i_gene][0] = sample_hits
        if enrichment_fold <= 1.0:
            continue

            
        # Get hypergeometric p-value
        p_value = hyper.getUpperCumulativeHyperP(
            int(total_size),
            int(sample_size),
            int(total_hits),
            int(sample_hits))
        
        p_values[i_gene][0] = p_value
        
        
    return p_values, enrichment_folds, module_hits



def calculate_enrichments_final(network, global_network, feature_matrix, HRR_sizes, global_hrr, 
                                recovered_terms, gene_ids, genes_to_predict, min_hits=2, only_global=False):
    
    gene_ind = dict(zip(gene_ids, list(range(len(gene_ids)))))
    total_size = network.shape[1]
    

    # number of genes for each GO
    feature_sizes = np.sum(feature_matrix, axis=0)                      # per col: sum over all rows
    

    # Prepare emtpy matrices for p-values and enrichment folds
    p_values = np.ones([len(genes_to_predict), feature_matrix.shape[1]], dtype=np.float32)
    enrichment_folds = np.ones([len(genes_to_predict), feature_matrix.shape[1]], dtype=np.float32)
    module_hits = np.ones([len(genes_to_predict), feature_matrix.shape[1]], dtype=np.float32)


    hyper = Hyper()

    # Loop genes
    for i_gene, gene in enumerate(genes_to_predict):
        
        # translate gene id into index (of common_ids, not the original compendium index!).
        ind = gene_ind[gene]                            

        
        for i_feature, feature_size in enumerate(feature_sizes):

            hrr_size = HRR_sizes[i_feature]

            if i_feature not in recovered_terms:
                
                if hrr_size != global_hrr:
                    raise ValueError('[ERROR] {}, {}, {}'.format(i_feature, hrr_size, global_hrr))

                assert hrr_size == global_hrr
                

                current_cluster = global_network[ind]
            elif only_global:
                continue
            else:
                # copy gene cluster 
                current_cluster = np.copy(network[ind])
                
                # cut-off cluster
                cutoff_network(current_cluster, hrr_size)
            
        
            # Get indexes of cluster genes
            cluster = np.logical_not(np.isnan(current_cluster))                # logical_not(isnan()): element-wise NaN -> False

            
            sample_size = np.sum(cluster)
            if sample_size == 0:
                raise ValueError('[ERROR] Empty cluster')
            elif sample_size == total_size:
                raise ValueError('[ERROR] cluster = all genes')

                
            # test if indices still match
            assert len(cluster) == len(feature_matrix)

            
            # number of cluster genes for this GO
            cluster_feature_hits = np.sum(feature_matrix[cluster, i_feature])
            assert isinstance(cluster_feature_hits, np.float64)


            # skip feature when too few hits
            if cluster_feature_hits < min_hits: 
                continue


            # Get enrichment fold
            total_hits = feature_size
            sample_hits = cluster_feature_hits

            enrichment_fold = (float(sample_hits)/sample_size) / (float(total_hits)/total_size)
            enrichment_folds[i_gene][i_feature] = enrichment_fold
            module_hits[i_gene][i_feature] = sample_hits
            if enrichment_fold <= 1.0:
                continue


            # Get hypergeometric p-value
            p_value = hyper.getUpperCumulativeHyperP(
                int(total_size),
                int(sample_size),
                int(total_hits),
                int(sample_hits))

            p_values[i_gene][i_feature] = p_value
        
        #break
        #if i_gene > 500:
        #    break


    return p_values, enrichment_folds, module_hits




def Fbeta(recall, precision, beta):
    return (1+beta**2)*recall*precision/(recall+(beta**2)*precision)



def evaluate_enrichments(p_values, enrichment_folds, gold_standard, ICu_matrix, go_ids, scanning=False, fdr=0.05, use_FC=False, use_Fmax=True, beta=1,
                            correct_MQSE=False, use_ICu=True, extend_after_predicting=False, reduce_final_prediction=True):
    '''
    Predicts functions with the specified parameters. Returns predictions and
    evaluation metrics.

    returns tuple:
        predictions
        recall
        precision
        F1/Fmax
        max allowed p-value
        predictions per gene
        highest IC of correct predictions per gene
        score_cutoff
    '''

    # Flatten the p-value matrix & sort
    p_values_flat = np.squeeze(p_values.reshape([1, -1]))
    p_values_flat[np.equal(p_values_flat, 0)] = 1e-25
    sorted_p_values = np.argsort(p_values_flat)


    # Find p-value cutoff at given FDR
    correction = np.sum(np.less(p_values_flat, 1))  # #pvalues < 1
    max_p = 0.

    for c, i in enumerate(sorted_p_values):
        
        if p_values_flat[i]*correction/(c+1) <= fdr:
            max_p = p_values_flat[i]
        else:
            break
    
    # Check if any predictions can be made:
    if max_p == 0:
        return np.zeros_like(p_values, dtype=int), 0., 0., 0., 0., 0., np.zeros([len(p_values)]), 0.

    # Generate a dummy ICu_matrix
    if use_ICu:
        dummy_ICu_matrix = ICu_matrix
    else:
        dummy_ICu_matrix = np.ones_like(ICu_matrix)


    # Get genes with features
    has_features = 1*np.greater(np.sum(gold_standard, axis=1, keepdims=True), 0)    # this is probably still equal to gold_standard (only when scanning HRR)
    if scanning:
        assert has_features.tolist() == gold_standard.tolist()



    # Check if Fmax is required
    if not use_Fmax:
        raise Exception('[EXCEPTION] not using Fmax has not been adapted to per GO optimization')


    # Prepare to find Fmax
    gold_standard_flat = np.squeeze((gold_standard*dummy_ICu_matrix).reshape([1, -1]))                 # this also has no effect as there is only 1 ICu value
    full_standard_flat = np.squeeze((np.ones_like(gold_standard)*dummy_ICu_matrix).reshape([1, -1]))   # this tags all genes??? what is this used for?
    has_features_flat = np.squeeze((np.ones_like(gold_standard)*has_features).reshape([1, -1]))
    total_positives = np.sum(gold_standard*dummy_ICu_matrix)                                           # this is (#current feature genes)*(current feature ICu) 
    total_predictions = 0.
    IC_predictions = 0.
    true_positives = 0.
    Fmax = 0.
    score_cutoff = 0.
    
    
    # Sort the predictions
    if use_FC:
        prediction_scores = np.squeeze(-enrichment_folds.reshape([1, -1]))
        sorted_predictions = np.argsort(prediction_scores)
    elif correct_MQSE:
        enrichment_folds_flat = np.squeeze(enrichment_folds.reshape([1, -1]))
        prediction_scores = np.power(10, enrichment_folds_flat*np.log10(p_values_flat))       # =LN
        sorted_predictions = np.argsort(prediction_scores)
        enrichment_folds_flat = None
    else:
        prediction_scores = p_values_flat
        sorted_predictions = sorted_p_values


    # Prepare predictions as a flat array of integers.
    predictions_flat = np.empty_like(p_values_flat, dtype=np.int)
    predictions_flat.fill(len(predictions_flat))


    # Loop through predictions (i.e. MQSE scores), updating F1 and finding Fmax
    total_positives_at_fmax = 0
    for i, prediction_index in enumerate(sorted_predictions):           # sorted_predictions: argsorted MQSE scores (thus: sorted genes)
        
        # Check if FDR is met
        if p_values_flat[prediction_index] >= max_p:                    # once this is the case, it's impossible to go further, right? -> break ipv continue
            continue

        # Set the current integer in the predictions_flat array
        predictions_flat[prediction_index] = total_predictions

        # Count the positive
        total_predictions += 1    # TP and FN !

        # Check if this gene had the current feature   -> only TP are allowed through. So here, this step should be disregarded as we still need FPs
        if not scanning:
            if not has_features_flat[prediction_index]:
                continue

        # Count IC of new prediction
        IC_predictions += full_standard_flat[prediction_index]      # this is TP + FP

        # Count if true positive
        true_positives += gold_standard_flat[prediction_index] # Will be zero if no true positives     # counts are ICu weighted 
        
        if true_positives == 0:
            continue
                
        # Check if next prediction has the same score
        next_score = 100.
        if i+1 < len(sorted_predictions):                                # avoid out of bound error
            next_score = prediction_scores[sorted_predictions[i+1]]

        # Only evaluate Fmax if the next prediction has a different score.
        # Otherwise, we may exclude some predictions with an identical score
        # to the last prediction made.
        if next_score != prediction_scores[prediction_index]:

            # Calcualte F1
            recall = float(true_positives)/total_positives
            precision = float(true_positives)/IC_predictions
            F1 = Fbeta(recall, precision, beta)
            
            if F1 > Fmax:
                Fmax = F1
                score_cutoff = prediction_scores[prediction_index]          # Ok so for this HRR, this MQSE score cutoff has best F1
                total_positives_at_fmax = total_predictions


    # At this point, 'total_positives_at_fmax tells us at which iteration we found
    # Fmax. This means that any integer in predictions_flat smaller or equal to this            # smaller or equal to? -> np.less ??
    # value indicates that this predictions is made at Fmax
    predictions = 1*(np.less(predictions_flat, total_positives_at_fmax)).reshape(p_values.shape)


    if extend_after_predicting:
        extend_with_go(predictions, go_ids, smaller_is_better=False)    # this is skipped


    # Now we simply recalculate evaluation metrics for these predictions
    true_positives = np.sum(predictions*gold_standard*dummy_ICu_matrix)
    if scanning:
        total_predictions = np.sum(predictions*dummy_ICu_matrix)                         # without has_features, not only TPs are counted (also FP). To correctly compare with a regular run, I should read in the OG has_features
    else:
        total_predictions = np.sum(has_features*predictions*dummy_ICu_matrix)
    total_positives = np.sum(gold_standard*dummy_ICu_matrix)

    if total_predictions == 0:
        return np.zeros_like(p_values, dtype=int), 0., 0., 0., 0., 0., np.zeros([len(p_values)]), 0.

    recall = float(true_positives)/total_positives
    precision = float(true_positives)/total_predictions
    F1 = Fbeta(recall, precision, beta)
    highest_ic_per_gene = np.max(predictions*gold_standard*ICu_matrix, axis=1)

    if reduce_final_prediction:
        reduce_with_go(predictions, go_ids)

    predictions_per_gene = np.sum(predictions, axis=1)

    return predictions, recall, precision, F1, max_p, predictions_per_gene, highest_ic_per_gene, score_cutoff
