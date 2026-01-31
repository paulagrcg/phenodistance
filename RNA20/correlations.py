import numpy as np 
from scipy.stats import pearsonr, spearmanr
import re
from pathlib import Path
import matplotlib.pyplot as plt

r_values = []
p_values = []
structures = []
h_mean_e_val_correlations = []  # Store r and p for H_mean vs e_val
h_vs_rho_by_sample = {}
structure_by_sample = {}
significant_samples = set()  # Track samples with significant H vs rho correlation

rna20_dir = "/home/pg520/phenodistance/RNA20"
rna20_files = sorted(Path(rna20_dir).glob("RNA20_site_scan_summary_*.txt"))

# First pass: identify significant samples
for rna_file in rna20_files:
    with open(rna_file, 'r') as f:
        content = f.read()
    
    # Extract number from filename
    match = re.search(r'RNA20_site_scan_summary_(\d+)\.txt', rna_file.name)
    if match:
        sample_num = match.group(1)
    
    # Extract correlation values (H vs rho)
    corr_match = re.search(r'r=([\d.e+-]+),\s*p=([\d.e+-]+)', content)
    if corr_match:
        r = float(corr_match.group(1))
        p = float(corr_match.group(2))
        
        if p < 0.05:
            significant_samples.add(sample_num)
            r_values.append(r)
            p_values.append(p)
            
            # Extract MFE Structure
            mfe_match = re.search(r'MFE Structure: (.+)', content)
            mfe = mfe_match.group(1) if mfe_match else "N/A"
            structures.append(mfe)
            h_vs_rho_by_sample[sample_num] = {
                'r': r,
                'p': p
            }
            structure_by_sample[sample_num] = mfe
            print(f"Significant Sample {sample_num}: r={r}, p={p}")

# Second pass: calculate H_mean vs e_val only for significant samples
for rna_file in rna20_files:
    with open(rna_file, 'r') as f:
        content = f.read()
    
    # Extract number from filename
    match = re.search(r'RNA20_site_scan_summary_(\d+)\.txt', rna_file.name)
    if match:
        sample_num = match.group(1)
    
    # Only process if this sample is significant in H vs rho
    if sample_num not in significant_samples:
        continue
    
    # Calculate correlation between H_mean and e_val
    lines = content.split('\n')
    h_mean_data = []
    e_val_data = []
    
    for line in lines:
        # Skip header and empty lines
        if line.startswith('#') or not line.strip():
            continue
        
        try:
            parts = line.split()
            if len(parts) >= 8:
                # Column 1 is H_mean, Column 7 is e_val
                h_mean = float(parts[1])
                e_val = float(parts[7])
                h_mean_data.append(h_mean)
                e_val_data.append(e_val)
        except (ValueError, IndexError):
            continue
    
    # Calculate correlation if we have data
    if len(h_mean_data) > 1:
        h_mean_corr, h_mean_p = pearsonr(h_mean_data, e_val_data)
        h_mean_e_val_correlations.append({
            'sample': sample_num,
            'r': h_mean_corr,
            'p': h_mean_p,
            'n_sites': len(h_mean_data)
        })

print(f"\nTotal samples with p < 0.05 (H vs rho): {len(r_values)}")
print(f"{'='*70}")
print("Correlation between H_mean and e_val (only for significant H vs rho samples):")
print(f"{'='*70}")

# Track samples with significant H_mean vs e_val correlation (p < 0.05)
significant_h_mean_e_val = []
for corr_data in sorted(h_mean_e_val_correlations, key=lambda x: int(x['sample'])):
    print(f"Sample {corr_data['sample']:>3}: r={corr_data['r']:>8.6f}, p={corr_data['p']:>12.6e}, n_sites={corr_data['n_sites']}")
    if corr_data['p'] < 0.05:
        significant_h_mean_e_val.append(corr_data['sample'])

print(f"\nSignificant H_mean vs e_val correlations (p < 0.05): {significant_h_mean_e_val}")

# Write CSV for significant H vs rho structures
h_vs_rho_csv = '/home/pg520/phenodistance/RNA20_structures_h_vs_rho_significant.csv'
with open(h_vs_rho_csv, 'w') as table_file:
    table_file.write('sample,structure,r,p\n')
    for sample in sorted(significant_samples, key=lambda x: int(x)):
        structure = structure_by_sample.get(sample, 'N/A')
        h_vs_rho = h_vs_rho_by_sample.get(sample, {})
        h_vs_rho_r = h_vs_rho.get('r', 'N/A')
        h_vs_rho_p = h_vs_rho.get('p', 'N/A')
        table_file.write(f"{sample},{structure},{h_vs_rho_r},{h_vs_rho_p}\n")

# Write CSV for significant H_mean vs e_val structures
h_mean_e_val_csv = '/home/pg520/phenodistance/RNA20_structures_h_vs_e_val_significant.csv'
h_mean_e_val_map = {item['sample']: item for item in h_mean_e_val_correlations}
with open(h_mean_e_val_csv, 'w') as table_file:
    table_file.write('sample,structure,r,p,n_sites\n')
    for sample in sorted(significant_h_mean_e_val, key=lambda x: int(x)):
        structure = structure_by_sample.get(sample, 'N/A')
        h_mean_vs_e_val = h_mean_e_val_map.get(sample, {})
        h_mean_vs_e_val_r = h_mean_vs_e_val.get('r', 'N/A')
        h_mean_vs_e_val_p = h_mean_vs_e_val.get('p', 'N/A')
        h_mean_vs_e_val_n = h_mean_vs_e_val.get('n_sites', 'N/A')
        table_file.write(f"{sample},{structure},{h_mean_vs_e_val_r},{h_mean_vs_e_val_p},{h_mean_vs_e_val_n}\n")

print(f"\nStructure table saved to {h_vs_rho_csv}")
print(f"Structure table saved to {h_mean_e_val_csv}")

# Plot histogram of Pearson r values (overlayed)
plt.figure(figsize=(10, 6))

h_vs_rho_r = r_values
h_mean_vs_e_val_r = [item['r'] for item in h_mean_e_val_correlations]

bins = 15

if h_mean_vs_e_val_r:
    plt.hist(
        h_mean_vs_e_val_r,
        bins=bins,
        density=True,
        color='red',
        edgecolor='black',
        alpha=0.6,
        label=r'$\langle H_p^{[i]}\rangle,\ e_p^{[i]}$'
    )

if h_vs_rho_r:
    plt.hist(
        h_vs_rho_r,
        bins=bins,
        density=True,
        color='blue',
        edgecolor='black',
        alpha=0.6,
        label=r'$\langle H_p^{[i]}\rangle,\ \rho_p^{[i]}$'
    )

plt.xlabel('Pearson r', fontsize=12)
plt.ylabel('Normalised Frequency', fontsize=12)
plt.title('Pearson r Distributions', fontsize=14)
plt.legend(fontsize=10)
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('/home/pg520/phenodistance/r_values_histogram.pdf', dpi=300)
print("\nHistogram saved to r_values_histogram.pdf")
plt.show()