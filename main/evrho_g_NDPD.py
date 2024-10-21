from functions.structurefunctions import *
import pickle 
import sys 
import os


if __name__ == "__main__":

        a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
        gpmap  = pickle.load(a_file)
        K = 4
        L = 12
        i = sys.argv[1]

        evg = defaultdict(float)
        rhog = defaultdict(float)

        f = open('/rds/user/pg520/hpc-work/seqsfiles/seqsfiles/file'+i+'.txt',"rb")
        lines = f.readlines()
        for line in lines:
                seq = line[0:L]
                seq = seq.decode('utf-8')
                evg[seq] = evolvabilitygND_PD(gpmap,seq,K,L)[seq]
                rhog[seq] = robustnessgND_PD(gpmap,seq,K,L)[seq]
        a_file = open("/rds/user/pg520/hpc-work/evgNDPD/evgNDPD"+i+".pkl","wb")
        pickle.dump(evg,a_file)
        if not os.path.exists('/rds/user/pg520/hpc-work/rhogNDPD'):
                os.makedirs('/rds/user/pg520/hpc-work/rhogNDPD')
        a_file = open("/rds/user/pg520/hpc-work/rhogNDPD/rhoNDPD"+i+".pkl","wb")
        pickle.dump(rhog,a_file)
