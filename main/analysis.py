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
    hamming_local_site, edgeD, phenos_site= hamming_local_D_PD_site(DGPmap)
    with open("/home/pg520/phenodistance/data/hamming_local_D_PD_site.pkl","wb") as f:
        pickle.dump(hamming_local_site,f)
    with open("/home/pg520/phenodistance/data/edge_hamming_local_D_PD_site.pkl","wb") as f:
        pickle.dump(edgeD,f)
    with open("/home/pg520/phenodistance/data/phenos_site.pkl","wb") as f:
        pickle.dump(phenos_site,f)
    hamming_local, edge, phenos, phenosevolvability, phenosevweighted = hamming_local_D_PD_site_nodel(DGPmap)
    with open("/home/pg520/phenodistance/data/hamming_local_site_nodel.pkl","wb") as f:
        pickle.dump(hamming_local,f)
    with open("/home/pg520/phenodistance/data/edgeD_nodel.pkl","wb") as f:
        pickle.dump(edge,f)
    with open("/home/pg520/phenodistance/data/phenos_site_nodel.pkl","wb") as f:
        pickle.dump(phenos,f)
    with open("/home/pg520/phenodistance/data/phenos_site_evol.pkl","wb") as f:
        pickle.dump(phenosevolvability,f)
    with open("/home/pg520/phenodistance/data/phenos_site_evol_weight.pkl","wb") as f:
        pickle.dump(phenosevweighted,f)
#    hammingglobalDPD_nodel = hamming_global_D_PD_nodel(neutralsets,L)