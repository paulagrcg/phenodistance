import itertools as it
from collections import Counter
import pickle
from collections import defaultdict
import time
import sys
import os
import math
import random 
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
                evgdictd[seq]+=stringdif(p,foldmut)*(1-pprime)
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
        phvsprobseq = extractnormalisedprobs(dictsuboptRNA12[seq],L)
        for phenotype,probg in phvsprobseq.values():
            foldList.append(phenotype)
            probsList.append(probg)
        resolutiongp[seq] = random.choices(foldList, weights = probsList, k=1)[0]
    return resolutiongp