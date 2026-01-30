import RNA
import numpy as np
from scipy.stats import pearsonr
from typing import Optional


def mutationalneighbours_site(seq, site):
    mutations = {'A': ['C','U','G'],'C': ['A','U','G'],'G': ['A','U','C'], 'U':['A','G','C']}
    return [seq[:site] + m + seq[site+1:] for m in mutations[str(seq[site])]]

def hamming(seq1, seq2):
    d = 0
    for i, j in zip(seq1, seq2):
        if i != j: 
            d += 1
    return d / len(seq1)

def mfe_structure(seq: str) -> str:
    fc = RNA.fold_compound(seq)
    ss, _ = fc.mfe()
    return ss

def hamming_str(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) / len(a)

def mutational_neighbours_site(seq: str, site: int):
    nts = ["A", "C", "G", "U"]
    wt = seq[site]
    for nt in nts:
        if nt != wt:
            yield seq[:site] + nt + seq[site+1:]

def site_scan_totals_non_neutral(sequence: str, wt_structure: Optional[str] = None):
    """
    For each site i:
      - H_total[i] = sum of Hamming distances WT->mutant over NON-neutral mutants only
      - n_non_neutral[i] = count of non-neutral mutants at that site
      - rho[i] = fraction neutral mutants at that site (out of 3)
      - e[i] = number of distinct phenotypes accessible at that site (incl WT if reached)
    """
    if wt_structure is None:
        wt_structure = mfe_structure(sequence)

    L = len(sequence)
    out = {}

    for i in range(L):
        H_total_list = []
        n_non_neutral = 0
        neutral = 0
        phenos_e = set()

        for mut in mutational_neighbours_site(sequence, i):
            ss = mfe_structure(mut)

            if ss == wt_structure:
                neutral += 1
            else:
                H_total_list.append(hamming_str(wt_structure, ss))
                n_non_neutral += 1
                phenos_e.add(ss)

        out[i] = {
            "H_total": np.sum(H_total_list),
            "H_total_list": H_total_list,
            "n_non_neutral": n_non_neutral,
            "rho": neutral / 3.0,
            "e phenos": phenos_e,
        }

    return out

def aggregate_RNA12_style_non_neutral(sequences: list[str], target_structure: Optional[str] = None):
    """
    RNA12-style phenotype-level aggregation using a finite sample:

    For each site i:
      ⟨H_{q[i]}⟩ = (sum over all sequences of H_total[i]) / (sum over all sequences of n_non_neutral[i])

    Neutral mutations are excluded from both sums.
    Sequences are optionally filtered to those with MFE == target_structure.
    """
    kept = []
    for seq in sequences:
        wt_ss = mfe_structure(seq)
        if target_structure is None or wt_ss == target_structure:
            kept.append((seq, wt_ss))

    if not kept:
        raise ValueError("No sequences matched the target structure (or empty list).")

    f_q = len(kept)
    L = len(kept[0][0])
    if any(len(s) != L for s, _ in kept):
        raise ValueError("All sequences must have the same length.")

    H_total_sum = np.zeros(L, dtype=float)
    n_non_neutral_sum = np.zeros(L, dtype=int)
    rho_vals = [[] for _ in range(L)]
    H_total = [[] for _ in range(L)]
    e_phenos_set = [set() for _ in range(L)]
    N_total_phenotypes = set()
    
    for seq, wt_ss in kept:
        stats = site_scan_totals_non_neutral(seq, wt_structure=wt_ss)
        for i in range(L):
            H_total_sum[i] += stats[i]["H_total"]
            H_total[i].extend(stats[i]["H_total_list"])
            n_non_neutral_sum[i] += stats[i]["n_non_neutral"]
            rho_vals[i].append(stats[i]["rho"])
            e_phenos_set[i].update(stats[i]["e phenos"])
            N_total_phenotypes.update(e_phenos_set[i])
    
    summary = {}
    for i in range(L):
        H_total_array = np.array(H_total[i])
        H_mean = (H_total_sum[i] / n_non_neutral_sum[i]) if n_non_neutral_sum[i] > 0 else 0.0
        H_mean_from_array = np.mean(H_total_array)
        H_std_from_array = np.std(H_total_array)
        e_vals = len(e_phenos_set[i]) / len(N_total_phenotypes)
        summary[i] = {
            "H_mean": float(H_mean),
            "H_mean_from_array": float(H_mean_from_array),
            "H_std_from_array": float(H_std_from_array),
            "H_total": float(H_total_sum[i]),
            "n_non_neutral": int(n_non_neutral_sum[i]),
            "rho_mean": float(np.mean(rho_vals[i])),
            "e_val": float(e_vals),
            "f_q": int(f_q),
        }

    return summary, kept

def sitewise_correlations(summary):
    sites = sorted(summary.keys())

    H = np.array([summary[i]["H_mean"] for i in sites])
    rho = np.array([summary[i]["rho_mean"] for i in sites])

    r_Hrho, p_Hrho = pearsonr(H, rho)

    return {
        "H_vs_rho": (r_Hrho, p_Hrho),
        "H": H,
        "rho": rho,
    }


def scan_sites_mfe(seq, samplesize):
    """
    Sequentially scan through sites to find neutral mutations that preserve MFE structure.
    
    Parameters:
    seq: initial RNA sequence
    samplesize: number of sequences to find in the neutral space
    
    Returns:
    seqs: list of sequences (all with the same MFE structure)
    structures: list of MFE structures (all identical)
    """
    
    # Get MFE structure of wild-type
    wt_structure = mfe_structure(seq)
    
    seqs = [seq]
    structures = [wt_structure]
    current_seq = seq
    samplesizecount = 1  # Start count at 1 since we already have the wild-type
    
    seqlen = len(seq)
    site = 0
    
    while samplesizecount < samplesize:
        site_neutral = False
        mutationchoices = mutationalneighbours_site(current_seq, site)
        
        # Try random mutations at this site until finding a neutral one
        while not site_neutral and mutationchoices:
            m = np.random.randint(0, len(mutationchoices))
            seqmut = mutationchoices.pop(m)
            
            # Get MFE structure of mutant
            mut_structure = mfe_structure(seqmut)
            
            # If MFE structure is preserved, keep this mutation
            if wt_structure == mut_structure:
                seqs.append(seqmut)
                structures.append(mut_structure)
                current_seq = seqmut
                site_neutral = True
                samplesizecount += 1
                #print(f"Found neutral sequence {samplesizecount}/{samplesize} at site {site}")
                
                if samplesizecount == samplesize:
                    break
                
                site += 1
                if site == seqlen:
                    site = 0
        
        # If no neutral mutation found at this site, move to next site
        if not site_neutral:
            site += 1
            if site == seqlen:
                site = 0
        
        if samplesizecount == samplesize:
            break
    
    print(f"\nTotal neutral sequences found: {len(seqs)}")
    print(f"MFE Structure: {wt_structure}")
    return seqs, structures

# Main execution
if __name__ == "__main__":
    import sys
    import time
    
    # Read sequences from file
    with open('samples_sequences.txt', 'r') as f:
        sequences = [line.strip() for line in f if line.strip()]
    
    # Get sequence at index from command line argument
    idx = int(sys.argv[1])
    seq = sequences[idx]
    
    q = mfe_structure(seq)
    start_time = time.time()
    seqs, structures = scan_sites_mfe(seq, samplesize=1e6)
    summary, kept = aggregate_RNA12_style_non_neutral(seqs, q)
    dictcorr = sitewise_correlations(summary)["H_vs_rho"]
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")
    np.savetxt(f"RNA20_site_scan_summary_{idx}.txt",
               np.array([[i,
                          summary[i]["H_mean"],
                          summary[i]["H_mean_from_array"],
                          summary[i]["H_std_from_array"],
                          summary[i]["H_total"],
                          summary[i]["n_non_neutral"],
                          summary[i]["rho_mean"],
                          summary[i]["e_val"],
                          summary[i]["f_q"]] for i in sorted(summary.keys())]),
               fmt=['%d'] + ['%.6f'] * 8,
               header="site H_mean H_mean_from_array H_std_from_array H_total n_non_neutral rho_mean e_val f_q")
    print(f"Correlation H vs rho: r={dictcorr[0]}, p={dictcorr[1]}")