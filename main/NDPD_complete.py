import pickle
from collections import defaultdict

if __name__ == "__main__":

     # evgNDPDtot = defaultdict(float)
     # rhogNDPDtot = defaultdict(float)
     # for i in range(0,167):
     #      with open("/rds/user/pg520/hpc-work/evgNDPD/evgNDPD"+str(i)+".pkl","rb") as f:
     #           evgNDPD = pickle.load(f)
     #           evgNDPDtot = {**evgNDPDtot,**evgNDPD}
     #      with open("/rds/user/pg520/hpc-work/rhogNDPD/rhogNDPD"+str(i)+".pkl","rb") as f:
     #           rhogNDPD = pickle.load(f)
     #           rhogNDPDtot = {**rhogNDPDtot,**rhogNDPD}

     # with open("/home/pg520/phenodistance/data/evgNDPD.pkl","wb") as f:
     #      pickle.dump(evgNDPDtot,f)
     # with open("/home/pg520/phenodistance/data/rhogNDPD.pkl","wb") as f:
     #      pickle.dump(rhogNDPDtot,f)

     evpNDPDtot = defaultdict(float)
     rhopNDPDtot = defaultdict(float)

     for i in range(0,271):
          with open("/rds/user/pg520/hpc-work/evpNDPD/evpNDPD"+str(i)+".pkl","rb") as f:
               evpNDPD = pickle.load(f)
               evpNDPDtot = {**evpNDPDtot,**evpNDPD}
          with open("/rds/user/pg520/hpc-work/evpNDPD/rhopNDPD"+str(i)+".pkl","rb") as f:
               rhopNDPD = pickle.load(f)
               rhopNDPDtot = {**rhopNDPDtot,**rhopNDPD}

     with open("/home/pg520/phenodistance/data/evpNDPD.pkl","wb") as f:
          pickle.dump(evpNDPDtot,f)
     with open("/home/pg520/phenodistance/data/rhopNDPD.pkl","wb") as f:
          pickle.dump(rhopNDPDtot,f)