#!/usr/bin/python
"""
This script executes different network inference algorithms and returns the
networks as numpy arrays that can be passed to the evaluation functions.
"""
import numpy as np
import pandas as pd
from sys import stderr
from random import seed, shuffle
import gc

#####################
# GENERATE NETWORKS #
#####################

def simple_correlation(compendium):
    """
    Calculates regular PCC values between all genes in a compendium. Returns
    a fully connected network

    # Input:
    compendium = np.array( [ [ float ] ] ) -> genes x samples  expression values

    # output:
    network = np.array( [ [  float ]  ] )  -> genes x genes    correlation coefficients
    """
    # Calculate expression values minus the means and the variances
    minmeans = compendium - np.mean(compendium, axis=1, keepdims=True)
    variances = np.sqrt(np.sum(np.square(minmeans), axis=1, keepdims=True))
    
    # Calculate correlations between all genes (covariance / product of variances)
    network = np.dot(minmeans, minmeans.T) / np.dot(variances, variances.T)
    return network


def weighted_correlation(compendium, weights):
    """
    Calculates regular PCC values between all genes in a compendium. Returns
    a fully connected network

    # Input:
    compendium = np.array( [ [ float ] ] ) -> genes x samples  expression values

    # output:
    network = np.array( [ [  float ]  ] )  -> genes x genes    correlation coefficients
    """
    # Calculate expression values minus the means and the variances
    minmeans = compendium - np.average(compendium, weights=weights, axis=1).reshape((compendium.shape[0],1))
    variances = np.sqrt(np.sum( weights*minmeans/np.sum(weights)*minmeans , axis=1, keepdims=True))

    # Calculate correlations between all genes (covariance / product of variances)
    network = np.dot(weights * minmeans / np.sum(weights), minmeans.T) / np.dot(variances, variances.T)
    return network

def clr(compendium):
    
    """
    """
    
    geo_means = np.tile(stats.mstats.gmean(compendium, axis=1), (compendium.shape[1], 1)).T
    cent_log_ratio = np.log(compendium/geo_means)
    
    return cent_log_ratio

def getVar(array):
    
    minmeans = array - np.mean(array, axis=1, keepdims=True)
    variances = np.sum(np.square(minmeans), axis=1, keepdims=True)/array.shape[1]

    return variances


def pairwiseRho(compendium, verbose=True):
    """
    Calculates Rho values between all genes in a compendium. Returns
    a fully connected network

    # Input:
    compendium = np.array( [ [ float ] ] ) -> genes x samples  expression values

    # output:
    network = np.array( [ [  float ]  ] )  -> genes x genes    proportionality values
    """
    
    # init network
    n_genes = compendium.shape[0]
    n_samples = compendium.shape[1]
    network = np.zeros((n_genes, n_genes))

    # var(A_j)
    variances = getVar(compendium)
    
    # Covariance
    minmeans = compendium - np.mean(compendium, axis=1, keepdims=True)
    covariances = np.dot(minmeans, minmeans.T)/n_samples
    
    for i in range(n_genes):

        if verbose:
            print('{}/{}'.format(i+1, n_genes), end='\r')
            if i % 100 == 0:
                sys.stdout.flush()

        # dummy atlas from current row
        repeated_row = np.tile(compendium[i], (n_genes-i, 1))

        # var(A_i)
        row_variance = np.tile(variances[i], (n_genes-i, 1))
        
        # var(A_i - A_j)
        diff_variances = row_variance+variances[i:]-2*covariances[i,i:].reshape((n_genes-i, 1))

        # 1-var(A_i - A_j)/(var(A_i) + var(A_j))
        rho = 1-diff_variances/(row_variance+variances[i:])
    
        # substitute current row
        network[i,i:] = rho.flatten()
    
    network[np.equal(network, 0)] = -2
    network = np.maximum(network, network.T)
    
    assert network.min() >= -1
    assert network.max() <= 1.000001
    
    return network


#######################
# NETWORK CONVERSIONS #
#######################


def convert_hrr(network, weights_are_ranks=False, reorder=False, use_MR=False):
    """
    Converts a fully connected network to HRR.
    """
    BATCH_SIZE = 1000 # To reduce memory usage, only sort for 1000 genes at a time.
    # Copy the original network, set edges between the same gene to max/min depending on high_to_low.
    # This ensures that the gene's perfect correlation to themselves will not shift the ranks by one
    # for all the other genes.
    hrr_network = np.copy(network)
    if weights_are_ranks:
        hrr_network[np.identity(len(hrr_network), dtype=np.bool)] = np.max(hrr_network)
    else:
        # set diagonal to minimum value of pcc array
        hrr_network[np.identity(len(hrr_network), dtype=np.bool)] = np.min(hrr_network)
    # Convert to ranks
    if not weights_are_ranks or reorder:
        for batch_start in range(0, len(network), BATCH_SIZE):
            batch_size = min(BATCH_SIZE, len(network)-batch_start)
            batch_ranks = np.argsort(hrr_network[batch_start:batch_start+batch_size], axis=1)
            # Reverse the sorting
            if not weights_are_ranks:
                batch_ranks = np.fliplr(batch_ranks)
            hrr_network[batch_start:batch_start+batch_size] = np.argsort(batch_ranks, axis=1) + 1 # +1 to make the lowest rank 1
    
    # Get the minimum ranks
    if use_MR:
        mr_network = (hrr_network * hrr_network.T)**(1/2)       # geometric mean of ranks
        return mr_network
    else:
        hrr_network = np.minimum(hrr_network, hrr_network.T)    # highest rank = lowest value
        return hrr_network


def convert_knn(network, cutoff=100, weights_are_ranks=False):
    """
    Converts a fully connected network to KNN. This network is not symmetrical!
    """
    BATCH_SIZE = 1000 # To reduce memory usage, only sort for 1000 genes at a time.

    # Copy the original network, set edges between the same gene to max/min depending on high_to_low.
    # This ensures that the gene's perfect correlation to themselves will not shift the ranks by one
    # for all the other genes.
    knn_network = network[:]
    if weights_are_ranks:
        knn_network[np.identity(len(knn_network), dtype=np.bool)] = np.max(knn_network)
    else:
        knn_network[np.identity(len(knn_network), dtype=np.bool)] = np.min(knn_network)

    # Loop batches
    for batch_start in range(0, len(network), BATCH_SIZE):
        batch_size = min(BATCH_SIZE, len(network)-batch_start)
        batch_ranks = np.argsort(knn_network[batch_start:batch_start+batch_size], axis=1)
        if not weights_are_ranks: # Reverse the sorting
            batch_ranks = np.fliplr(batch_ranks)
        knn_network[batch_start:batch_start+batch_size] = np.argsort(batch_ranks, axis=1) + 1 # +1 to make the lowest rank 1

    cutoff_network(knn_network, threshold=cutoff, weights_are_ranks=True)

    return knn_network



def cutoff_network(network, threshold=100, weights_are_ranks=True):
    """
    Replaces values above the threshold with 'nan'.
    In case of correlations as weights, values below the threshold will be removed instead.
    """

    if weights_are_ranks:
        network[np.isnan(network)] = threshold + 1 # Get rid of numpy warning
        network[np.greater(network, threshold)] = float('nan')
    
    else:
        network[np.isnan(network)] = threshold - 1 # Get rid of numpy warning
        network[np.less(network, threshold)] = float('nan')



###################
# READ/WRITE DATA #
###################

# Compendium

def read_compendium(compendium_file, expressed_only=True, cutoff=2):
    """
    Input compendia files are tab-separated and formatted as follows:

    #      sample1  sample2  sample3  ...
    gene1    1.2      3.0      2.8    ...
    gene2    9.4      4.6      2.7    ...
    gene3    1.1      7.6      3.6    ...
     ...     ...      ...      ...    ...

    This function reads in such a file and returns the sample ids, gene ids as lists
    and the expression values as a numpy array

    # Input:
    compendium_file = "file_path"

    # Output:
    sample_ids = [ sample_id  ]
    gene_ids   = [  gene_id   ]
    compendium = np.array([ [ values ] ]) -> genes x samples
    """
    sample_ids = None
    gene_ids = []
    gene_set = set()
    compendium = []
    with open(compendium_file, 'r') as reader:
        sample_ids = reader.readline().strip().split("\t")[1:]
        line_counter = 1
        for line in reader:
            line_counter += 1
            split = line.strip().split("\t")
            if len(split) != len(sample_ids)+1:
                stderr.write("[ERROR] compendium file badly formatted at line "+str(line_counter)+"\n")
                stderr.write("[ERROR]     -> Should be "+str(len(sample_ids)+1)+" columns long\n")
                exit()
            gene_id = split[0]
            if gene_id in gene_set:
                stderr.write("[ERROR] compendium file contains a duplicated entry at line "+str(line_counter)+"\n")
                stderr.write("[ERROR]     -> Identifier of "+str(gene_id)+" appears twice\n")
                exit()
            try:
                values = list(map(float, split[1:]))
                if max(values) == min(values):              # skip if no variation
                    continue
                if expressed_only and max(values) < cutoff: # if requested, skip unexpressed
                    continue
                gene_ids.append(gene_id)
                compendium.append(values)
            except ValueError:
                stderr.write("[ERROR] compendium file badly formatted at line "+str(line_counter)+"\n")
                stderr.write("[ERROR]     -> Non-numerical value in expression data\n")
                exit()
    compendium = np.array(compendium, dtype=np.float32)
    return compendium, gene_ids, sample_ids

def write_compendium(compendium, gene_ids, sample_ids, output_file):
    with open(output_file, 'w') as writer:
        writer.write("#\t"+"\t".join(sample_ids)+"\n")
        for i_gene in range(len(gene_is)):
            writer.write(gene_ids[i_gene]+"\t"+"\t".join(map(str, compendium[i_gene]))+"\n")

def read_compendium_npz(compendium_file):
    container = np.load(compendium_file)
    return container['compendium'], container['gene_ids'], container['sample_ids']

def write_compendium_npz(output_file, compendium, gene_ids, sample_ids):
    np.savez(output_file, compendium=compendium, gene_ids=gene_ids, sample_ids=sample_ids)

def shuffle_compendium_genes(compendium, gene_ids):
    """
    Shuffles genes in the compendium to remove the sorted gene bias.
    """
    seed(1990)
    shuf_gene_ind = list(range(len(gene_ids)))
    shuffle(shuf_gene_ind)
    shuf_gene_ids = [gene_ids[i] for i in shuf_gene_ind]
    shuf_compendium = compendium[shuf_gene_ind]
    return shuf_compendium, shuf_gene_ids

def write_network_tsv(input_file, hrr_cutoff, output_file, gene_ids=None):

    # read network    
    if not gene_ids:
        hrr_network, gene_ids = read_network_npz(input_file)
        assert hrr_network.shape == (len(gene_ids), len(gene_ids))
    else:
        hrr_network = input_file
    gene_ids = np.array(gene_ids)
    
    # cutoff network    
    hrr_network[np.greater(hrr_network, hrr_cutoff)] = float('nan')
    hrr_clusters = np.logical_not(np.isnan(hrr_network))

    # fetch edges
    seeds, prey, hrr_values = [], [], []
    for (i, row), bool_row in zip(enumerate(hrr_network), hrr_clusters):
        
        current_cluster_genes = gene_ids[bool_row]
        current_hrr_values = row[bool_row]
        
        seeds += [gene_ids[i]]*len(current_cluster_genes)
        prey += list(current_cluster_genes)
        hrr_values += [int(el) for el in list(current_hrr_values)]
    
    # store
    hrr_network_long = pd.DataFrame({'seed': seeds, 'prey': prey, 'HRR': hrr_values})
    hrr_network_long = hrr_network_long.sort_values(['seed','HRR','prey'])
    hrr_network_long.to_csv(output_file, sep='\t', index=False)


def read_network_npz(network_file):
    container = np.load(network_file, encoding='bytes')
    return container['network'], container['gene_ids']

def write_network_npz(output_file, network, gene_ids):
    np.savez(output_file, network=network, gene_ids=gene_ids)
