from collections import Counter
from collections import defaultdict
import random 
import sys
import functools
sys.path.append('/home/pg520/phenodistance/main/functions/')
from basefunctions import *
import numpy as np
import os 
import pickle

stn = {'A': '0', 'U': '1', 'C': '2', 'G': '3'}

def evodyn(initpop, T, phenoq, L=12, gpmap = {}, mu = 0.05, generations = 10000, Npop = 100):
        popgenos = {}
        #initpop is a tupple with initial target pheno and initial population list
        popgenos = np.empty((generations+1, Npop),dtype=int)
        timediscovered = defaultdict(int)
        hammingdist = defaultdict(int)
        fitness = defaultdict(float)
        phenos = defaultdict(list)

        for i in range(1,generations+1): # i is generation number
                if i ==1:
                        genos = np.array([initpop[1]]*Npop)
                        fitness[phenoq] = 1
                else: genos = newpop
                popgenos[i] = [int(''.join([stn[n] for n in geno])) for geno in genos] #taken from i-1 mutated population
                unique, counts = np.unique(popgenos[i], return_counts=True)
                freq = dict(zip(unique, counts)) #counts the number of times each genotype appears in the population
                phenos = []
                for j in range(0,Npop):
                        phenos[i].append(gpmap[genos[j]])

                for pheno in phenos:
                        if pheno not in list(fitness.keys()):
                                fitness[pheno] = 1+np.heaviside(i - T, hamming(pheno,phenoq)) #fitness increases by Hamming distance once generation T is reached
         
        
                #mean fitness of population at gen i       
                meanfitness[i] = sum(list(popgen.values()))/len(popgen)
                #data collection of target probabilities
                if i == generations:
                    probstarg0end = probstarg0
                    probstarg1end = probstarg1
                #probstarg0dict[i] = probstarg0
                #probstarg1dict[i] = probstarg1
                meanprobs0[i-1] = np.mean(list(probstarg0.values()))
                meanprobs1[i-1] = np.mean(list(probstarg1.values()))
                meanprobs0seqscom[i-1] =np.mean(probstarg0seqscom) 
                meanprobs1seqscom[i-1] =np.mean(probstarg1seqscom)
                #roulette wheel selection
                newpop = []
                choices = []
                index = []
                fsum = sum(list(popgen.values()))
                probRWS = defaultdict(float)

                for g,ph in popgen.items():
                        probRWS[g] = float(ph/fsum)
                        choices.append(g[0])

                #sampling of new population 
                newpop = [random.choices(choices, weights=list(probRWS.values()), k=1)[0] for num in range(0,Npop)]
                     
                #mutations for population i+1
                newpop = [mutation(seq, mu) for seq in newpop]

        #return meanfitness,seqstargetlist,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0dict,probstarg1dict,meanprobs0,meanprobs1
        return meanfitness,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0end,probstarg1end,meanprobs0seqscom,meanprobs1seqscom,meanprobs0,meanprobs1,alphalist
if __name__ == "__main__":
        pair = int(sys.argv[1])
        fitlands = int(sys.argv[2])#1.random 2.hamming 3.random with inverse
        samplenum = int(sys.argv[3])
        plasticoption = int(sys.argv[4])
        gengap = float(sys.argv[5])
        generations = int(sys.argv[6])
        mu = float(sys.argv[7])
        minprob = int(sys.argv[8])
        minprob = minprob/100.
        maxprob = int(sys.argv[9])
        maxprob = maxprob/100.
        initoption = int(sys.argv[10]) #extremeoption -> 0: soft extremes, 1: true extremes, 2: random from seqscommon
        extremeoption = int(sys.argv[11]) #e.g. if we want to start simulation with target 0, the initial population should be at extreme of target 1
        Npop = int(sys.argv[12])
        run = int(sys.argv[13])
        samples = 100 
        path = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/"    
        pathx = "/rds/user/pg520/hpc-work/targetflipping/pair_"+str(pair)+"/regimedata/"
        path1 = pathx + "fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+ "/"
        if not os.path.isdir(pathx): os.mkdir(pathx)
        if not os.path.isdir(path): os.mkdir(path)
        if not os.path.isdir(path1): os.mkdir(path1)

        with open("/home/pg520/target_flipping/phenopairs.pkl","rb") as f: phenopairs = pickle.load(f)
        phenopair = phenopairs[pair]
        pheno0 = phenopair[0]
        pheno1 = phenopair[1]
        del phenopairs 
        

        with open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl","rb") as f: gpmap = pickle.load(f)

        seqscommon = np.load(path + "seqscommon.npy")
        #probs_common_0 = np.load(path + "probs_common_0.npy")
        #probs_common_1 = np.load(path + "probs_common_1.npy")

        #seqstarget = np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy")
        seq_soft_extreme_0 = np.load(path + "seqssoftextreme0.npy")
        seq_soft_extreme_1 = np.load(path + "seqssoftextreme1.npy")
        seqstrueextreme0 = stringtoint(np.load(path + "seqstrueextreme0.npy"))
        seqstrueextreme1 = stringtoint(np.load(path + "seqstrueextreme1.npy"))
        seqscommon = stringtoint(np.load(path + "seqscommon.npy"))
        #seqstarget = stringtoint(np.load(path + "seqstarget_"+str(minprob)+"_"+str(maxprob)+".npy"))
        
        #seq_extreme_0_max = np.load(path + "seq_trueextreme_0_max.npy")
        #seq_extreme_1_max = np.load(path + "seq_trueextreme_1_max.npy")
        #probs_extreme_0_max = np.load(path + "/probs_trueextreme_0_max.npy")
        #probs_extreme_1_max = np.load(path + "/probs_trueextreme_1_max.npy")


        #with open("/rds/user/pg520/hpc-work/entropy.pkl","rb") as f: entropy = pickle.load(f)
        #with open("/rds/user/pg520/hpc-work/evgND.pkl","rb") as f: evgND = pickle.load(f)
        #with open("/rds/user/pg520/hpc-work/rhogND.pkl","rb") as f: rhogND = pickle.load(f)
        entropy = defaultdict(float)
        evgND = defaultdict(float)
        rhogND = defaultdict(float)
        
        #initial population
        #print("initoption:", initoption)
        #print("extremeoption:", extremeoption)

        initpop = []
        if initoption == 0:
                if extremeoption == 0: initpop = (1,seq_soft_extreme_0)
                if extremeoption == 1: initpop = (0,seq_soft_extreme_1)
        elif initoption == 1: 
                if extremeoption == 0: initpop = (1,seq_extreme_0_max)
        elif initoption == 1: 
                if extremeoption == 0: initpop = (1,seq_extreme_0_max)
                if extremeoption == 1: initpop = (0,seq_extreme_1_max)
        elif initoption == 2:
                p0val = 1
                p1val = 0
                while val == True:
                        seq = random.choice(seqscommon)
                        indexseq = list(seqscommon).index(seq)
                        if seq not in seqstarget: val == False
                initpop = (random.choice([0,1]),seq)
        #print('targetoption and genotype start:', initpop)
        #meanrobust, meanentropy, meanfitness, adaptedphenonum, targetgenosgen, probstarg0dict, probstarg1dict, listseqprobs, popgenos, meanevolgnd = evodyn(rhogND=rhogND, evgND=evgND, entropy=entropy,targets=phenopair,initpop=initpop, genostargetset=seqstarget, seqscommon= seqscommon, samplenum = samplenum, plasticoption = plasticoption, gengap = gengap, L=12, gpmap = gpmap, mu = mu, generations = generations, Npop = Npop, fitlands = fitlands)
        dfs = []
        column_mean = defaultdict(list)
        column_std = defaultdict(list)
        probstot0 = []
        probstot1 = []
        for sample in range(0,samples):
                meanfitness,seqscommonlist,extreme0list,extreme1list,otherlist,probstarg0end,probstarg1end,meanprobs0seqscom,meanprobs1seqscom,meanprobs0,meanprobs1,alphalist = evodyn(rhogND=rhogND, evgND=evgND, entropy=entropy,targets=phenopair,initpop=initpop, seqscommon= seqscommon, extreme0 = seqstrueextreme0, extreme1 = seqstrueextreme1, samplenum = samplenum, plasticoption = plasticoption, gengap = gengap, L=12, gpmap = gpmap, mu = mu, generations = generations, Npop = Npop, fitlands = fitlands)
                dictsampleprobs ={'meanfitness': list(meanfitness.values()),'seqscommon': seqscommonlist,'extreme0': extreme0list, 'extreme1': extreme1list, 'other':otherlist,'meanprobs0seqscom':meanprobs0seqscom, 'meanprobs1seqscom':meanprobs1seqscom, 'meanprobs0':meanprobs0, 'meanprobs1':meanprobs1,'alphalist': alphalist}
                #print(dictsampleprobs)
                #dictsample = {'meanfitness': list(meanfitness.values()), 'seqscommon': seqscommonlist, 'seqstarget': seqstargetlist, 'extreme0': extreme0list, 'extreme1': extreme1list, 'other':otherlist}
                #df = pd.DataFrame(dictsample)
                #dfs.append(df)
                df = pd.DataFrame(dictsampleprobs)
                dfs.append(df)
                for v0,v1 in zip(list(probstarg0end.values()), list(probstarg1end.values())):
                    probstot0.append(v0)
                    probstot1.append(v1)
                
    
        concatenated_df = pd.concat(dfs)
        dfmean = concatenated_df.groupby(concatenated_df.index).mean()
        dfstd = concatenated_df.groupby(concatenated_df.index).std()
        np.save(path1 + "probstarg0.npy", np.array(probstot0))
        np.save(path1 + "probstarg1.npy", np.array(probstot1))
        #dfmean.to_csv(path1 + 'meandata'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
        #dfstd.to_csv(path1 + 'stddata'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
        dfmean.to_csv(path1 + 'meanprobsalphadata_'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
        dfstd.to_csv(path1 + 'stdprobsalphadata_'+"fitlands"+str(fitlands)+"_initoption"+str(initoption)+"_extremeoption"+str(extremeoption)+"_samplenum"+str(samplenum)+"_plasticoption"+str(plasticoption)+"_gengap"+str(gengap)+"_mu"+str(mu)+"_generations"+str(generations)+"_Npop"+str(Npop)+'.csv')
        #with open(path1 + "probstarg0dict.pkl","wb") as f:
        #        pickle.dump(probstarg0dict, f)
        #with open(path1 + "probstarg1dict.pkl","wb") as f:
        #        pickle.dump(probstarg1dict, f)