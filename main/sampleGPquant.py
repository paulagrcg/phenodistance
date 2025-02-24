import pickle
import os
from collections import defaultdict
from collections import Counter
from functions.structurefunctions import neutralsets_DPD

def get_files_in_directory(directory):
    files = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            files.append(filename)
    return files


if __name__ == "__main__":
    L = 12
    K = 4
    directory_path = "/rds/user/pg520/hpc-work/sampleGP/"
    files = get_files_in_directory(directory_path)
    
    results = defaultdict(dict)
    for file in files:
        file_path = os.path.join(directory_path, file)
        a_file = open(file_path, "rb")
        gpmap  = pickle.load(a_file)
        folds = list(gpmap.values())
        neutralsets = Counter(folds)
        result = neutralsets_DPD(neutralsets,L)
        results[file] = result

    averageresult = defaultdict(float)
    print(len(results))
    for result in results:
        for fold,nss in result.items():
            averageresult[fold] += nss/len(results)
    with open("averagenss_NDPD.pkl", "wb") as f:
        pickle.dump(averageresult, f)

