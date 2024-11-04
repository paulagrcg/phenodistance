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
 
#    hamminglocalDPD_nodel, edgeD = hamming_local_D_PD_nodel(DGPmap)
#    hammingglobalDPD_nodel = hamming_global_D_PD_nodel(neutralsets,L)
#    phipq = phipqD(DGPmap,neutralsets,K,L)
#    with open("/home/pg520/phenodistance/data/hamminglocalDPD_nodel.pkl","wb") as f:
#        pickle.dump(hamminglocalDPD_nodel,f)
#    with open("/home/pg520/phenodistance/data/edgeD.pkl","wb") as f:
#        pickle.dump(edgeD,f)
#    with open("/home/pg520/phenodistance/data/hammingglobalDPD_nodel.pkl","wb") as f:
#        pickle.dump(hammingglobalDPD_nodel,f)
#    with open("/home/pg520/phenodistance/data/phipq.pkl","wb") as f:
#        pickle.dump(phipq,f)
    hamming_local_site, edgeD = hamming_local_D_PD_site(DGPmap)
    with open("/home/pg520/phenodistance/data/hamming_local_D_PD_site.pkl","wb") as f:
        pickle.dump(hamming_local_site,f)
    hamming_global_site = hamming_global_D_PD_site(DGPmap, neutralsets)
    with open("/home/pg520/phenodistance/data/hamming_global_D_PD_site.pkl","wb") as f:
        pickle.dump(hamming_global_site,f)
#    hammingglobalDPD_nodel = hamming_global_D_PD_nodel(neutralsets,L)