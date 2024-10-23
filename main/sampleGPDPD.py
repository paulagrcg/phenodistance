from collections import Counter
from functions.structurefunctions import *
import pickle

if __name__ == "__main__":
    K = 4
    L = 12
    i = sys.argv[1]

    a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
    gpmap  = pickle.load(a_file)
    sampleGPmap = sampleGP(gpmap,L)

    neutralsets = Counter(list(sampleGPmap.values()))

    rho_g_pd,rho_p_pd = robustnessD_PD(sampleGPmap,folddict,K,L)

    with open("/home/pg520/phenodistance/data/rhopDPD.pkl","wb") as f:
        pickle.dump(rho_p_pd,f)

    with open("/home/pg520/phenodistance/data/rhogDPD.pkl","wb") as f:
        pickle.dump(rho_g_pd,f)

    ev_g_pd = evolvabilitygD_PD(sampleGPmap)
    ev_p_pd = evolvabilitypD_PD(sampleGPmap)

    with open("/home/pg520/phenodistance/data/evgDPD.pkl","wb") as f:
        pickle.dump(ev_g_pd,f)
    with open("/home/pg520/phenodistance/data/evpDPD.pkl","wb") as f:
        pickle.dump(ev_p_pd,f)

    neutralsetsDPD = neutralsets_DPD(neutralsets,L)

    with open("/home/pg520/phenodistance/data/neutralsetsDPD.pkl","wb") as f:
        pickle.dump(neutralsetsDPD,f)