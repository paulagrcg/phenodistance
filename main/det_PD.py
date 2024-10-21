import pickle
from collections import Counter
from functions.structurefunctions import *

if __name__ == "__main__":
    K = 4
    L = 12

    a_file = open("/home/pg520/phenodistance/data/DGPmap.pkl", "rb")
    DGPmap  = pickle.load(a_file)

    with open("/rds/user/pg520/hpc-work/folddictt.pkl", "rb") as f:
        folddict = pickle.load(f)

    rho_g_pd,rho_p_pd = robustnessD_PD(DGPmap,folddict,K,L)

    with open("/home/pg520/phenodistance/data/rhopDPD.pkl","wb") as f:
        pickle.dump(rho_p_pd,f)

    with open("/home/pg520/phenodistance/data/rhogDPD.pkl","wb") as f:
        pickle.dump(rho_g_pd,f)

    ev_g_pd = evolvabilitygD_PD(DGPmap)
    ev_p_pd = evolvabilitypD_PD(DGPmap)

    with open("/home/pg520/phenodistance/data/evgDPD.pkl","wb") as f:
        pickle.dump(ev_g_pd,f)
    with open("/home/pg520/phenodistance/data/evpDPD.pkl","wb") as f:
        pickle.dump(ev_p_pd,f)