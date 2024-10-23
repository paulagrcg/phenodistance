import pickle
from collections import Counter
from functions.structurefunctions import *
import os 

def calculate_det_quant(gpmap,K,L):
    folds = list(gpmap.values())
    neutralsets = Counter(folds)
    rho_g_pd,rho_p_pd = robustnessD_PD(gpmap,folds,K,L)

    ev_g_pd = evolvabilitygD_PD(gpmap)
    ev_p_pd = evolvabilitypD_PD(gpmap)

    neutralsetsDPD = neutralsets_DPD(neutralsets,L)

    return rho_g_pd,rho_p_pd,ev_g_pd,ev_p_pd,neutralsetsDPD

if __name__ == "__main__":
    K = 4
    L = 12
    i = int(sys.argv[1]) #if i = 0 then we do MFE 
    
    if i > 0:
        a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
        gpmap  = pickle.load(a_file)

        sampleGPmap = sampleGP(gpmap,L)
        rho_g_pd,rho_p_pd,ev_g_pd,ev_p_pd,neutralsetsDPD = calculate_det_quant(sampleGPmap,K,L)
        
        if not os.path.exists('/rds/user/pg520/hpc-work/samples'):
            os.makedirs('/rds/user/pg520/hpc-work/samples')
        with open("/rds/user/pg520/hpc-work/samples/rhopDPD"+str(i)+".pkl","wb") as f:
            pickle.dump(rho_p_pd,f)
        with open("/rds/user/pg520/hpc-work/samples/rhogDPD"+str(i)+".pkl","wb") as f:
            pickle.dump(rho_g_pd,f)
        with open("/rds/user/pg520/hpc-work/samples/evgDPD"+str(i)+".pkl","wb") as f:
            pickle.dump(ev_g_pd,f)
        with open("/rds/user/pg520/hpc-work/samples/evpDPD"+str(i)+".pkl","wb") as f:
            pickle.dump(ev_p_pd,f)
        with open("/rds/user/pg520/hpc-work/samples/neutralsetsDPD"+str(i)+".pkl","wb") as f:
            pickle.dump(neutralsetsDPD,f)
    else: 
        #MFE gp map
        a_file = open("/home/pg520/phenodistance/data/DGPmap.pkl", "rb")
        gpmap  = pickle.load(a_file)
        rho_g_pd,rho_p_pd,ev_g_pd,ev_p_pd,neutralsetsDPD = calculate_det_quant(gpmap,K,L)
        with open("/home/pg520/phenodistance/data/rhopDPD.pkl","wb") as f:
            pickle.dump(rho_p_pd,f)
        with open("/home/pg520/phenodistance/data/rhogDPD.pkl","wb") as f:
            pickle.dump(rho_g_pd,f)
        with open("/home/pg520/phenodistance/data/evgDPD.pkl","wb") as f:
            pickle.dump(ev_g_pd,f)
        with open("/home/pg520/phenodistance/data/evpDPD.pkl","wb") as f:
            pickle.dump(ev_p_pd,f)
        with open("/home/pg520/phenodistance/data/neutralsetsDPD.pkl","wb") as f:
            pickle.dump(neutralsetsDPD,f)

    