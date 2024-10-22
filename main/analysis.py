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
 
    hamminglocalDPD = hamming_local_D_PD(DGPmap)
    hammingglobalDPD = hamming_global_D_PD(neutralsets,L)
    phipq = phipqD(DGPmap,neutralsets,K,L)
    with open("/home/pg520/phenodistance/data/hamminglocalDPD.pkl","wb") as f:
        pickle.dump(hamminglocalDPD,f)
    with open("/home/pg520/phenodistance/data/hammingglobalDPD.pkl","wb") as f:
        pickle.dump(hammingglobalDPD,f)
    with open("/home/pg520/phenodistance/data/phipq.pkl","wb") as f:
        pickle.dump(phipq,f)
