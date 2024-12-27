import numpy as np 
import subprocess
import os
import math
from functions.basefunctions import *
from collections import defaultdict
from functionalRNA import *
import sys
import time
import pickle
from scipy.stats import gaussian_kde


def plasticitydistance(fRNAprob1, fRNAprob2):
    distances = {}
    for x, y,key1,key2 in zip(fRNAprob1.values(), fRNAprob2.values(), fRNAprob1.keys(), fRNAprob2.keys()):
        if key1 == key2:
            distance = np.sqrt((x - 0.5)**2 + (y - 0.5)**2)
            distances[key1] = distance
    return distances

def hamming_plasticity_optimal(fRNAhammingdistance, fRNAfolds, distances):
    scores = {}
    for name in fRNAhammingdistance.keys(): 
        scores[name] = fRNAhammingdistance[name] - distances[name] 
    # Remove entries from scores where the sequence of dots is in fRNAfolds
    filtered_scores = {key: value for key, value in scores.items() if '.' * len(key[1]) not in fRNAfolds[key]}

    # Select the data points with the highest scores
    selected_data = sorted(filtered_scores.items(), key=lambda item: item[1], reverse=True)

    # Extract the Hamming values and their corresponding minimum distances
    selected_hamming_values = [fRNAhammingdistance[name] for name, score in selected_data]
    selected_min_distances = [distances[name] for name, score in selected_data]
    selected_names = [name for name, score in selected_data]

    return selected_hamming_values, selected_min_distances, selected_names

def mutatesite(seq, site):
    mutations = {'A': ['C', 'U', 'G'], 'C': ['A', 'U', 'G'], 'G': ['A', 'U', 'C'], 'U': ['A', 'G', 'C']}
    return [seq[:site] + m + seq[site+1:] for m in mutations[str(seq[site])]]

def scan_sites(seq, samplesize):
    seqs = []
    folds1 = []
    folds2 = []
    probs1 = []
    probs2 = []
    samplesizecount = 0
    samplesizebreak = False

    seqoutput = suboptfolding(seq)
    seqs.append(seq)

    fold1 = seqoutput[0][0]
    fold2 = seqoutput[1][0]
    prob2 = float(seqoutput[1][2])
    prob1 = float(seqoutput[0][2])

    folds1.append(fold1)
    folds2.append(fold2)
    probs1.append(prob1)
    probs2.append(prob2)

    seqlen = len(seq)
    while samplesizecount < samplesize:
        for site in range(seqlen):
            for seqmut in mutatesite(seq, site):
                mutoutput = suboptfolding(seqmut)
                if fold1 == mutoutput[0][0] and fold2 == mutoutput[1][0]:  # if MFE is the same (so in the same neutral space)
                    #print(f"Match found: {fold1} == {mutoutput[0][0]} and {fold2} == {mutoutput[1][0]}")
                    seqs.append(seqmut)
                    fold1 = mutoutput[0][0]
                    fold2 = mutoutput[1][0]
                    prob2 = float(mutoutput[1][2])
                    prob1 = float(mutoutput[0][2])

                    folds1.append(fold1)
                    folds2.append(fold2)
                    probs1.append(prob1)
                    probs2.append(prob2)

                    samplesizecount += 1
                    seq = seqmut
                    samplesizebreak = True
                    break
            if samplesizebreak:
                #print(f"Breaking inner loop at site {site}")
                samplesizebreak = False
                break
        if samplesizebreak:
            #print("Breaking outer loop")
            break

    print(f"Total sequences found: {samplesizecount}")
    return seqs, folds1, folds2, probs1, probs2

def compute_pvalue(p1, p2, kde, data):
        # Evaluate the density at the point (p1, p2)
        point_density = kde([p1, p2])

        # Compute densities for all points
        all_densities = kde(data)

        # Calculate the p-value as the fraction of points with density <= point_density
        p_value = np.sum(all_densities <= point_density) / len(all_densities)

        return p_value

if __name__ == "__main__":

    
    with open('./data/fRNAhammingdistance.pkl','rb') as f:
        fRNAhammingdistance = pickle.load(f)
    with open('./data/fRNAprob2.pkl','rb') as f:
        fRNAprob2 =  pickle.load(f)
    with open('./data/fRNAprob1.pkl','rb') as f:
        fRNAprob1 =  pickle.load(f)
    with open('./data/fRNAfolds.pkl','rb') as f:
        fRNAfolds =  pickle.load(f)

    distances = plasticitydistance(fRNAprob1, fRNAprob2)
    selected_hamming_values, selected_min_distances, selected_names = hamming_plasticity_optimal(fRNAhammingdistance, fRNAfolds, distances)

    #site-scanning
    samplesize = int(sys.argv[1])
    seqposition = int(sys.argv[2])
    #i = 0
    #start = time.time()
    #seqs, folds1, folds2, probs1, probs2 = scan_sites(selected_sequences[i], samplesize)
    #end = time.time()

    #p1 and p2 for p-value calculation
    p1 = fRNAprob1[selected_names[seqposition]]
    p2 = fRNAprob2[selected_names[seqposition]]

    seqs, folds1, folds2, probs1, probs2 = scan_sites(selected_names[seqposition][1], samplesize)
    #p-value calculation
    data = np.vstack([probs1, probs2])
    kde = gaussian_kde(data)
    p_value = compute_pvalue(p1, p2, kde, data)

    #print(f"Time taken for site scanning: {end-start}")
    with open(f"./data/site_scanning_probs_pval_{seqposition}.pkl","wb") as f:
        pickle.dump({'probs1': probs1, 'probs2': probs2, 'p_value': p_value})
    
