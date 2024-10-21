from ..functions.structurefunctions import *

a_file = open("/rds/user/pg520/hpc-work/dictRNA12tot.pkl", "rb")
gpmap  = pickle.load(a_file)

i = sys.argv[1]
K = 4
L = 12

evg = defaultdict(float)

f = open('/rds/user/pg520/hpc-work/seqsfiles/seqsfiles/file'+i+'.txt',"rb")
lines = f.readlines()
for line in lines:
        seq = line[0:L]
        seq = seq.decode('utf-8')
        evg[seq] = evolvabilitygND_PD(gpmap,seq,4,12)[seq]
a_file = open("/rds/user/pg520/hpc-work/evgNDPD/evgNDPD"+i+".pkl","wb")
pickle.dump(evg,a_file)
