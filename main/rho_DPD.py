from ..functions.structurefunctions import *
import pickle

K = 4
L = 12

a_file = open("/home/pg520/phenodistance/data/DGPmap.pkl", "rb")
DGPmap  = pickle.load(a_file)

with open("/rds/user/pg520/hpc-work/folddictt.pkl", "rb") as f:
	folddict = pickle.load(f)

rho_g_pd,rho_p_pd = robustnessD_PD(DGPmap,folddict,K,L)

with open("/home/pg520/phenodistance/data/rhopDND.pkl","wb") as f:
    pickle.dump(rho_p_pd,f)

with open("/home/pg520/phenodistance/data/rhogDND.pkl","wb") as f:
    pickle.dump(rho_g_pd,f)