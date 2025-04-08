import pickle
K = 4
L = 12

#robustness and evolvability
from collections import Counter
from collections import defaultdict
from itertools import product
import random

def normalize(probs):
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]

def mutationalneighbours(seq):
    mutations = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    return [seq[:j] + m + seq[j+1:] for j in range(0,len(seq)) for m in mutations[str(seq[j])]]

def mutationalneighbours_site(seq,site):
    mutations = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    return [seq[:site] + m + seq[site+1:] for m in mutations[str(seq[site])]]

def extractnormalisedprobs(pboltzlist,L):
    probsnorm = []
    for p in pboltzlist:
        probsnorm.append(float(p[L+2:]))
    prob = normalize(probsnorm)
    return {pboltzlist[pi][0:L]: prob[pi] for pi in range(0,len(pboltzlist))}

def mutation(seq,mu):
    basemutation = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    for n in range(0,len(seq)):
        basemut = random.choices([1,0], weights=[mu,1-mu],k=1)[0]
        if basemut==1: #mutate a base
            mut = random.choices(basemutation[seq[n]],k=1)[0]
            p = seq[:n] + str(mut) + seq[n+1:]
            seq = p
    return seq


def hamming_vec(a,b):
        return sum(x!=y for x,y, in zip(list(a),list(b)))

def similarity(seq1,seq2):
    d = 0
    for i,j in zip(seq1,seq2):
        if i==j: 
            d+=1
    return d/len(seq1)

def hamming(seq1,seq2):
    d = 0
    for i,j in zip(seq1,seq2):
        if i!=j: 
            d+=1
    return d/len(seq1)
def robustnessD_PD(d,folds,K,L):
    s = list(d.keys())
    a = Counter(folds)
    rho_g_pd = {}
    for A in d.keys():
        rho_g_pd[A[0:L]] = 0
    rho_p_pd = {}
    for A in a.keys():
        rho_p_pd[A[0:L]] = 0
    for seq in d.keys():
        for mut in mutationalneighbours(seq):
                rho_g_pd[seq] += similarity(d[seq],d[mut])/((K-1)*L)

    for fold,count in a.items():
        for j in range(0,len(folds)):
            if fold[0:L]==folds[j][0:L]:
                rho_p_pd[fold[0:L]]+=rho_g_pd[s[j]]
        rho_p_pd[fold[0:L]]/= count

    return rho_g_pd,rho_p_pd

def evolvabilitygD_PD(gpmapdet):
    evgdictd = defaultdict(float)
    for seq,p in gpmapdet.items():
            prodfold = defaultdict(lambda:1)
            for newmutation in mutationalneighbours(seq):
                    foldmut=gpmapdet[newmutation]
                    if foldmut != p:
                        prodfold[foldmut] *=0
                    else: prodfold[foldmut]*=1
            for foldmut,pprime in prodfold.items():
                evgdictd[seq]+=hamming(p,foldmut)*(1-pprime)
    return evgdictd

def evolvabilitypD_PD(dictRNA12):
    evolp = defaultdict(float)
    prodfold = defaultdict(lambda:1) # I am assuming this is the product part - I have changed it so that this collects the data for all p&p' pairs
    for seq, phenotype in dictRNA12.items(): #by arranging the loops in this way, the slow part (extractnormalisedprobs(gpmap[seq],L)) is executed fewer times
        for newmutation in mutationalneighbours(seq):
            phenotypeprime = dictRNA12[newmutation]
            if phenotypeprime!= phenotype:
                prodfold[(phenotype, phenotypeprime)] *=0

    for (phenotype, phenotypeprime), val in prodfold.items():
        evolp[phenotype] +=hamming(phenotype,phenotypeprime)*(1-val)

    return evolp

def robustnessD(DGPmap,folds,K,L):   
    s = list(DGPmap.keys())
    a = Counter(folds)
    rho_g = {}
    for A in DGPmap.keys():
        rho_g[A[0:L]] = 0
    rho_p = {}
    for A in a.keys():
        rho_p[A[0:L]] = 0
    #genotypic robustness
    for seq in DGPmap.keys():
        for mut in mutationalneighbours(seq):
                if (DGPmap[mut]==DGPmap[seq]):rho_g[seq]+=1/((K-1)*L)
                else: continue
    #phenotypic robustness
    for fold,count in a.items():
        for j in range(0,len(folds)):
            if fold[0:L]==folds[j][0:L]:
                rho_p[fold[0:L]]+=rho_g[s[j]]
        rho_p[fold[0:L]]/= count
    
    return rho_g,rho_p

def evolvabilitygD(DGPmap):
    evgdictd = defaultdict(float)
    for seq,p in DGPmap.items():
            prodfold = defaultdict(lambda:1)
            for newmutation in mutationalneighbours(seq):            
                    foldmut=DGPmap[newmutation]
                    if foldmut != p:
                        prodfold[foldmut] *=0
                    else: prodfold[foldmut]*=1
            for pprime in prodfold.values():
                evgdictd[seq]+=(1-pprime)
    return evgdictd

def evolvabilitypD(DGPmap):
    evolp = defaultdict(float)
    prodfold = defaultdict(lambda:1) 
    for seq, phenotype in DGPmap.items():
        for newmutation in mutationalneighbours(seq):                    
            phenotypeprime = DGPmap[newmutation]
            if phenotypeprime!= phenotype:
                prodfold[(phenotype, phenotypeprime)] *=0

    for (phenotype, phenotypeprime), val in prodfold.items():
        evolp[phenotype] += (1-val)

    return evolp

with open("../data/DGPmap.pkl","rb") as f:
    DGPmap = pickle.load(f)


evgPD = evolvabilitygD_PD(DGPmap)
evpPD = evolvabilitypD_PD(DGPmap)
rho_g_pd, rho_p_pd = robustnessD_PD(DGPmap, DGPmap.keys(), K, L)

with open("../data/evgPD.pkl","wb") as f:
    pickle.dump(evgPD,f)
with open("../data/evpPD.pkl","wb") as f:
    pickle.dump(evpPD,f)
with open("../data/rho_g_pd.pkl","wb") as f:
    pickle.dump(rho_g_pd,f)
with open("../data/rho_p_pd.pkl","wb") as f:
    pickle.dump(rho_p_pd,f)