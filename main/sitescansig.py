import glob
import pandas as pd
import re
import numpy as np
from scipy.stats import norm
from collections import defaultdict
import pickle
import numpy as np

def find_diagonal_distances_val(prob1,prob2):
    distance = np.sqrt((prob1 - 0.5)**2 + (prob2 - 0.5)**2)/(np.sqrt(2)/2)
    return distance

data_folder = '../data/'

# Load the fRNAprob1 and fRNAprob2 dictionaries
with open(data_folder + 'fRNAprob1.pkl', 'rb') as f:
    fRNAprob1 = pickle.load(f)
with open(data_folder + 'fRNAprob2.pkl', 'rb') as f:
    fRNAprob2 = pickle.load(f)
# Load the fRNAhammingdistance dictionary
with open(data_folder + 'fRNAhammingdistance.pkl', 'rb') as f:
    fRNAhammingdistance = pickle.load(f)

# Load the names_to_analyse dictionary
data_folder = '../data/datatoanalyse/'
with open(data_folder + 'to_analyse.pkl', 'rb') as f:
    names_to_analyse = pickle.load(f)

# Pattern to match the files
file_pattern = data_folder + 'to_analyse*_ssize10000.pkl'

# List to store p-values from all files
fRNAsitescount = defaultdict(list)

files = glob.glob(file_pattern)
listtoanalyse = []
# Loop through each file and extract probs2
for file in files:
    match = re.search(r'to_analyse(\d+)_ssize10000', file)
    
    if match:
        seq_number = int(match.group(1))
        d_obs = find_diagonal_distances_val(fRNAprob1[tuple(names_to_analyse[seq_number])], fRNAprob2[tuple(names_to_analyse[seq_number])])

        # Read the file into a DataFrame
        with open(file, 'rb') as f:
            data = pickle.load(f)
        probs1 = list(data['probs1'])
        probs2 = list(data['probs2'])
        if (len(probs1), len(probs2)) == (10000, 10000):
            listtoanalyse.append(seq_number)
            dataset = np.array([find_diagonal_distances_val(p1, p2) for p1, p2 in zip(probs1,probs2)])
            mean = np.mean(dataset)
            std_dev = np.std(dataset)

        z_score = (d_obs - mean) / std_dev

        # One-tailed p-value (upper-tail)
        p_value_upper = 1 - norm.cdf(z_score)
        fRNAsitescount[(seq_number,tuple(names_to_analyse[seq_number]))] = [mean, std_dev,d_obs, p_value_upper]

with open('./data/fRNAsitesfile.pkl', 'wb') as f:
    pickle.dump(fRNAsitescount, f)