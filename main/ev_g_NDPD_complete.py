import pickle
from collections import defaultdict

evgNDPDtot = defaultdict(float)

for i in range(0,167):
    with open("/rds/user/pg520/hpc-work/evgNDPD/evgNDPD"+str(i)+".pkl","rb") as f:
         evgNDPD = pickle.load(f)
         evgNDPDtot = {**evgNDPDtot,**evgNDPD}

with open("/home/pg520/phenodistance/data/evgNDPD.pkl","wb") as f:
     pickle.dump(evgNDPDtot,f)

evpNDPDtot = defaultdict(float)

for i in range(0,271):
    with open("/rds/user/pg520/hpc-work/evpNDPD/evpNDPD"+str(i)+".pkl","rb") as f:
         evpNDPD = pickle.load(f)
         evpNDPDtot = {**evpNDPDtot,**evpNDPD}

with open("/home/pg520/phenodistance/data/evpNDPD.pkl","wb") as f:
     pickle.dump(evpNDPDtot,f)