import pickle
from collections import Counter
from functions.structurefunctions import *

if __name__ == "__main__":
    K = 4
    L = 12

    a_file = open("/home/pg520/phenodistance/data/DGPmap.pkl", "rb")
    DGPmap  = pickle.load(a_file)
    a_file = open("/home/pg520/phenodistance/data/neutralsets.pkl", "rb")
    neutralsets = pickle.load(a_file)
    
    phi_pq = phipqD(DGPmap,neutralsets,K,L)
    phi_pq_site = phipqD_site(DGPmap,neutralsets,K,L)

    with open('/home/pg520/phenodistance/data/phi_pq.pkl', 'wb') as f:
        pickle.dump(phi_pq, f)
    with open('/home/pg520/phenodistance/data/phi_pq_site.pkl', 'wb') as f:
        pickle.dump(phi_pq_site, f)
        
    """ 
    hamming_local_mean, hamming_local_std, edgediffnodel, phenos, phenosevolvability, phenosevweighted, robustnesssite, totrobust, edgenondel = hamming_local_D_PD_site_nodel(DGPmap)
    siteshammingmean = defaultdict(functools.partial(defaultdict, float))
    siteshammingstd = defaultdict(functools.partial(defaultdict, float))
    sitesrobustness = defaultdict(functools.partial(defaultdict, float))
    sitesevolvability = defaultdict(functools.partial(defaultdict, float))
    datarobust = defaultdict(functools.partial(defaultdict, float))
    for pheno,hammingmean in hamming_local_mean.items():
        for site, mean in hammingmean.items():
            siteshammingmean[pheno][site] = mean
            siteshammingstd[pheno][site]  = hamming_local_std[pheno][site]

            sitesrobustness[pheno][site] = robustnesssite[pheno][site]/edgenondel[pheno][site] #mutations that lead to p over total mutatations except del
            sitesevolvability[pheno][site] = phenosevolvability[pheno][site] / (len(neutralsets) - 1)

    with open('/home/pg520/phenodistance/data/sitesevolvability.pkl', 'wb') as f:
        pickle.dump(sitesevolvability, f)
    with open('/home/pg520/phenodistance/data/sitesrobustness.pkl', 'wb') as f:
        pickle.dump(sitesrobustness, f)
    with open('/home/pg520/phenodistance/data/siteshammingmean.pkl', 'wb') as f:
        pickle.dump(siteshammingmean, f)
    with open('/home/pg520/phenodistance/data/siteshammingstd.pkl', 'wb') as f:
        pickle.dump(siteshammingstd, f)
    """