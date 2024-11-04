from collections import Counter
from collections import defaultdict
import random 
import sys
import functools
sys.path.append('/home/pg520/phenodistance/main/functions/')
from basefunctions import *

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

def robustnessgND_PD(gpmap,seq,K,L):
         rhogndict = defaultdict(float)
         #first normalize probabilites:
         phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
         neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
         for phenotype, probg in phvsprobseq.items():
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut, probpgmut in phvsprobmut.items():
                    rhogndict[seq] += similarity(phenomut,phenotype)*probpgmut*probg/((K-1)*L) 
         del neighbourvsphvsprob

         return rhogndict
     
def NDfolds_setsizes(gpmap,K,L):
    
    folddict = defaultdict(list)

    for seq in gpmap.keys():
        phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
        for phenotype,probg in phvsprobseq.items():
            #folds.append(phenotype)
            #NDsetsize[phenotype] += probg
            folddict[phenotype].append([seq,probg])
            
    return folddict 


def robustnesspND_PD(folddict,gpmap):
    rhopndict = defaultdict(float)
    
    for f,seqprobs in folddict.items():
        probgtot = 0
        for seq_p in seqprobs:
            seq=seq_p[0]
            probg =float(seq_p[1])
            probgtot += float(seq_p[1])
            neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
            
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probpgmut in phvsprobmut.items():
                    #if f == phenomut:
                        rhopndict[f] += similarity(f,phenomut)*probpgmut*probg
            del neighbourvsphvsprob

        rhopndict[f]/=probgtot
        
    return rhopndict
    

def robustnesspND_PD0(folddictkey,folddictval,gpmap):
        f = folddictkey
        seqprobs = folddictval
        probgtot = 0
        rhop = 0
        for seq_p in seqprobs:
            seq=seq_p[0]
            probg =float(seq_p[1])
            probgtot += float(seq_p[1])
            neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}

            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probpgmut in phvsprobmut.items():
                    #if f == phenomut:
                        rhop += similarity(f,phenomut)*probpgmut*probg
            del neighbourvsphvsprob

        rhop/=probgtot

        return rhop


def evolvabilitygND_PD(gpmap,seq,K,L):
         evgndict = defaultdict(float)
         phvsprobseq = extractnormalisedprobs(gpmap[seq],L)
         neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
         for phenotype, probg in phvsprobseq.items():
            probfold = defaultdict(lambda:1)
            evgndictp = 0
            for newmutation, phvsprobmut in neighbourvsphvsprob.items():
                for phenomut,probpmut in phvsprobmut.items():
                    if phenotype != phenomut:
                       probfold[phenomut]*=(1-probpmut)  #probfold[(phenomut,newmutation)] *=(1-probpmut)
                    else: continue
            for t, prob in probfold.items():
                evgndictp+=hamming(phenotype,t)*(1-prob)
            evgndict[seq]+=evgndictp*probg
            
         del neighbourvsphvsprob
     
         return evgndict

def evolvabilitypND_PD0(folddictkey,folddictval,gpmap): 

    evolp = defaultdict(float)

    f = folddictkey
    probgprime = defaultdict(lambda:1)
    for seq_p in folddictval:
        seq=seq_p[0]
        probg =float(seq_p[1])
        neighbourvsphvsprob = {newmutation: extractnormalisedprobs(gpmap[newmutation],L) for newmutation in mutationalneighbours(seq)}
        
        for newmutation, phvsprobmut in neighbourvsphvsprob.items():
            for phenomut,probpgmut in phvsprobmut.items():
                if f != phenomut:
                    probgprime[phenomut] *=(1-probpgmut*probg)
                else: continue
        del neighbourvsphvsprob

    for phenomut,prob in probgprime.items():
        evolp[f]+=hamming(f,phenomut)*(1-prob)
           
    return evolp
    
def sampleGP(dictsuboptRNA12, L):
    resolutiongp = {}
    for seq,subopt in dictsuboptRNA12.items():
        foldList = []
        probsList = []
        phvsprobseq = extractnormalisedprobs(subopt,L)
        for phenotype,probg in phvsprobseq.items():
            foldList.append(phenotype)
            probsList.append(probg)
        resolutiongp[seq] = random.choices(foldList, weights = probsList, k=1)[0]
    return resolutiongp

def neutralsets_DPD(neutralsets,L):
    fPD= defaultdict(float)
    for key1 in neutralsets.keys():
        p1 = key1[:12]  
        for key2 in neutralsets.keys():
            p2 = key2[:12]  
            for d in range(L + 1):
                if similarity(p1, p2) == d / L:
                    fPD[p1] += neutralsets[key2] * (d / L)
    return fPD

def hamming_local_D_PD(gpmap):
    edge = defaultdict(float)
    hamming_local = defaultdict(functools.partial(defaultdict, float))
    for seq in gpmap.keys():
        for mut in mutationalneighbours(seq):
            if gpmap[seq] != gpmap[mut]: #ignore robustness term
                edge[gpmap[seq]] +=1
                hamming_local[gpmap[seq]][hamming(gpmap[seq],gpmap[mut])] +=1
    return hamming_local, edge

def hamming_local_D_PD_site(gpmap):
    edge = defaultdict(float)
    hamming_local = defaultdict(functools.partial(defaultdict, float))
    for seq in gpmap.keys():
        for site in range(0,len(seq)):
            for mut in mutationalneighbours_site(seq, site):
                if gpmap[seq] != gpmap[mut]: #ignore robustness term
                    edge[gpmap[seq]][site] +=1
                    hamming_local[gpmap[seq]][site] += hamming(gpmap[seq],gpmap[mut])
    return hamming_local, edge


def hamming_local_D_PD_nodel(gpmap):
    edge = defaultdict(float)
    L=12
    hamming_local = defaultdict(functools.partial(defaultdict, float))
    for seq in gpmap.keys():
        if gpmap[seq] == '.'*L: continue
        for mut in mutationalneighbours(seq):
            if gpmap[mut] == '.'*L or gpmap[seq] == gpmap[mut]: continue
            edge[gpmap[seq]] +=1
            hamming_local[gpmap[seq]][hamming(gpmap[seq],gpmap[mut])] +=1
    return hamming_local, edge

def hamming_global_D_PD(neutralsets,L):
    hamming_global = defaultdict(functools.partial(defaultdict, float))
    for fold in neutralsets.keys():
        fold = fold[0:L]
        for fold1 in neutralsets.keys():
            if fold1[0:L] != fold: #ignore robustness term
                foldd = fold1[0:L]
                hamming_global[fold][hamming(fold,foldd)] += neutralsets[fold1]
    return hamming_global

def hamming_global_D_PD_site(gpmap, neutralsets):
    edge = defaultdict(float)
    hamming_global = defaultdict(functools.partial(defaultdict, float))
    for seq in gpmap.keys():
        for site in range(0,len(seq)):
            for mut in mutationalneighbours_site(seq, site):
                if gpmap[seq] != gpmap[mut]: #ignore robustness term
                    hamming_global[gpmap[seq]][site] += hamming(gpmap[seq],gpmap[mut])/neutralsets[gpmap[mut]]
    return hamming_global

def hamming_global_D_PD_nodel(neutralsets,L):
    hamming_global = defaultdict(functools.partial(defaultdict, float))
    for fold in neutralsets.keys():
        fold = fold[0:L]
        if fold=='.'*L: continue
        for fold1 in neutralsets.keys():
            if fold1[0:L]=='.'*L or fold1[0:L] == fold: continue
            foldd = fold1[0:L]
            hamming_global[fold][hamming(fold,foldd)] += neutralsets[fold1]
    return hamming_global

def phipqD(gpmap,neutralsets,K,L):
    phi_pq = defaultdict(functools.partial(defaultdict, float))
    for seq in gpmap.keys():
        if seq == '.'*L: continue
        for mut in mutationalneighbours(seq):
            phi_pq[gpmap[seq]][gpmap[mut]] +=1/(neutralsets[gpmap[seq]+'\n']*(K-1)*L)
    return phi_pq
