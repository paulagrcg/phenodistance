import pickle
from collections import defaultdict
import sys
import glob
import os

def average_dicts(files):
    avg = defaultdict(float)
    nfiles = len(files)
    if nfiles == 0:
        print("No files found.")
        return avg
    for fname in files:
        with open(fname, 'rb') as handle:
            data = pickle.load(handle)
            for k, v in data.items():
                avg[k] += v
    for k in avg:
        avg[k] /= nfiles
    return avg

if __name__ == "__main__":

    i = int(sys.argv[1])

    if i == 0:
        files = sorted(glob.glob("/rds/user/pg520/hpc-work/samples/neutralsetsDPD*.pkl"))
        print(f"Processing {len(files)} files for neutral sets")
        avnss = average_dicts(files)
        with open("/home/pg520/phenodistance/data/avnssDPD.pkl", "wb") as f:
            pickle.dump(avnss, f)

    if i == 1:
        files = sorted(glob.glob("/rds/user/pg520/hpc-work/samples/rhogDPD*.pkl"))
        print(f"Processing {len(files)} files for rho_g")
        avrhog = average_dicts(files)
        with open("/home/pg520/phenodistance/data/avrhogDPD.pkl", "wb") as f:
            pickle.dump(avrhog, f)

    if i == 2:
        files = sorted(glob.glob("/rds/user/pg520/hpc-work/samples/rhopDPD*.pkl"))
        print(f"Processing {len(files)} files for rho_p")
        avrhop = average_dicts(files)
        with open("/home/pg520/phenodistance/data/avrhopDPD.pkl", "wb") as f:
            pickle.dump(avrhop, f)

    if i == 3:
        files = sorted(glob.glob("/rds/user/pg520/hpc-work/samples/evolgDPD*.pkl"))
        print(f"Processing {len(files)} files for evolution of g")
        avevolg = average_dicts(files)
        with open("/home/pg520/phenodistance/data/avevolgDPD.pkl", "wb") as f:
            pickle.dump(avevolg, f)

    if i == 4:
        files = sorted(glob.glob("/rds/user/pg520/hpc-work/samples/evolpDPD*.pkl"))
        print(f"Processing {len(files)} files for evolution of p")
        avevolp = average_dicts(files)
        with open("/home/pg520/phenodistance/data/avevolpDPD.pkl", "wb") as f:
            pickle.dump(avevolp, f)