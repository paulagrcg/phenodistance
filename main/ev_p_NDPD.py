from ..functions.structurefunctions import *

a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
gpmap  = pickle.load(a_file)

i = sys.argv[1]
K = 4
L = 12

with open("/rds/user/pg520/hpc-work/folddictt.pkl", "rb") as f:
	folddict = pickle.load(f)

evp = defaultdict(float)
folddictkeys = list(folddict.keys())[int(i)]
folddictvals = list(folddict.values())[int(i)]
evp = {folddictkeys: evolvabilitypND_PD0(folddictkeys,folddictvals,gpmap)}
a_file = open("/rds/user/pg520/hpc-work/evpNDPD/evpNDPD"+str(i)+".pkl","wb")
pickle.dump(evp,a_file)