import pickle
from collections import defaultdict
import sys 

if __name__ == "__main__":
  
    i = 1

    if i == 0:
        avnss = defaultdict(float)

        for j in range(1,101):
            with open("/rds/user/pg520/hpc-work/samples/neutralsetsDPD"+str(j)+".pkl", 'rb') as handle:
                nsets = pickle.load(handle)
                for k,nss in nsets.items():
                    avnss[k]+=nss/100.0

        with open("/home/pg520/phenodistance/data/avnssDPD.pkl","wb") as f:
            pickle.dump(avnss,f)

    if i == 1:
        avrhog = defaultdict(float)

        for j in range(1,251):
            with open("/rds/user/pg520/hpc-work/samples/rhogDPD"+str(j)+".pkl", 'rb') as handle:
                rho_gs = pickle.load(handle)
                for k,rhog in rho_gs.items():
                    avrhog[k]+=rhog/250.0

        with open("/home/pg520/phenodistance/data/avrhogDPD.pkl","wb") as f:
            pickle.dump(avrhog,f)
    if i == 2:
        avrhop = defaultdict(float)

        for j in range(1,251):
            with open("/rds/user/pg520/hpc-work/samples/rhopDPD"+str(j)+".pkl", 'rb') as handle:
                rho_ps = pickle.load(handle)
                for k,rhog in rho_ps.items():
                    avrhop[k]+=rhog/250.0

        with open("/home/pg520/phenodistance/data/avrhopDPD.pkl","wb") as f:
            pickle.dump(avrhop,f)

    if i == 3:
        avevolg = defaultdict(float)

        for j in range(1,251.0):
            with open("/rds/user/pg520/hpc-work/samples/evolgDPD"+str(j)+".pkl", 'rb') as handle:
                evol_gs = pickle.load(handle)
                for k,evolg in evol_gs.items():
                    avevolg[k]+=evolg/250.0

        with open("/home/pg520/phenodistance/data/avevolgDPD.pkl","wb") as f:
            pickle.dump(avevolg,f)

    if i == 4:
        avevolp= defaultdict(float)

        for j in range(1,251):
            with open("/rds/user/pg520/hpc-work/samples/evolpDPD"+str(j)+".pkl", 'rb') as handle:
                evol_ps = pickle.load(handle)
                for k,evolp in evol_ps.items():
                    avevolp[k]+=evolp/250.0

        with open("/home/pg520/phenodistance/data/avevolpDPD.pkl","wb") as f:
            pickle.dump(avevolp,f)
