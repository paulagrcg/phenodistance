import numpy as np 
import pickle 
from gpmapfunctions import *
from collections import defaultdict,Counter
import sys

#hamming_vec = np.vectorize(hamming_vec)
def hamming(a,b):
        return sum(x!=y for x,y, in zip(list(a),list(b)))/len(a)

def seqs_probs_total(folddict,pheno):
	seqsprobs =  np.array(folddict[pheno])
	seqs = list(seqsprobs[:,0])
	probs = list(np.float_(list(seqsprobs[:,1])))
	return seqs,probs

	
def localfreq(pheno, folddict):
	seqs, probs = seqs_probs_total(folddict,pheno)
	pheno_neighbours = defaultdict(float)
	del folddict[pheno]
	for pheno1 in (folddict.keys()):
		seqs1,probs1 = seqs_probs_total(folddict,pheno1)
		seqscount1 = Counter(seqs1)
		for seq in seqs:
			neighbours = mutationalneighbours(seq)
			for neigh in neighbours: 
				n = seqscount1[neigh]
				if n == 1:
					prob1 = probs1[np.where(np.array(seqs1) == neigh)[0][0]]
					pheno_neighbours[hamming(pheno,pheno1)] += prob1			
	return pheno_neighbours
if __name__ == "__main__":
	
	phenonum = int(sys.argv[1])	
	with open("/rds/user/pg520/hpc-work/folddictt.pkl","rb") as f: folddict = pickle.load(f)
	pheno = list(folddict.keys())[phenonum]

	localfreq = localfreq(pheno,folddict)
	with open("/rds/user/pg520/hpc-work/phenodistance/localfreq_"+str(phenonum)+".pkl", "wb") as f: pickle.dump(localfreq,f)

