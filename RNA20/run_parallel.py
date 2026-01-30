#!/usr/bin/env python3
"""
Run data.py for tasks 0-99 in parallel using multiprocessing.
"""
import sys
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
import RNA

# Import functions from data.py
from data import (
    mfe_structure, 
    scan_sites_mfe, 
    aggregate_RNA12_style_non_neutral,
    sitewise_correlations
)


def process_task(task_id, sequences):
    """Process a single task ID."""
    try:
        seq = sequences[task_id]
        
        # Check for invalid nucleotides
        valid_nts = set('ACGU')
        if not all(nt in valid_nts for nt in seq):
            invalid = set(seq) - valid_nts
            return {
                'task_id': task_id,
                'success': False,
                'error': f"Invalid nucleotides: {invalid}"
            }
        
        q = mfe_structure(seq)
        start_time = time.time()
        seqs, structures = scan_sites_mfe(seq, samplesize=1e5)
        summary, kept = aggregate_RNA12_style_non_neutral(seqs, q)
        dictcorr = sitewise_correlations(summary)["H_vs_rho"]
        end_time = time.time()
        
        # Save results
        np.savetxt(f"RNA20_site_scan_summary_{task_id}.txt",
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
        
        return {
            'task_id': task_id,
            'success': True,
            'time': end_time - start_time,
            'correlation': dictcorr
        }
    except Exception as e:
        return {
            'task_id': task_id,
            'success': False,
            'error': str(e)
        }


if __name__ == '__main__':
    import multiprocessing
    
    # Read sequences from file
    print("Loading sequences...")
    with open('sampled_sequences.txt', 'r') as f:
        sequences = [line.strip() for line in f if line.strip()]
    
    # Task IDs to process
    task_ids = range(0, 100)
    
    # Number of parallel workers (adjust based on your CPU)
    n_workers = multiprocessing.cpu_count()
    print(f"Running {len(task_ids)} tasks using {n_workers} workers...")
    
    # Run in parallel with progress bar
    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(process_task, task_id, sequences): task_id 
                   for task_id in task_ids}
        
        # Collect results with progress bar
        completed_futures = tqdm(as_completed(futures), total=len(task_ids)) if HAS_TQDM else as_completed(futures)
        for future in completed_futures:
            result = future.result()
            results.append(result)
            
            if result['success']:
                print(f"\nTask {result['task_id']} completed in {result['time']:.2f}s")
                print(f"  Correlation: r={result['correlation'][0]:.4f}, p={result['correlation'][1]:.4e}")
            else:
                print(f"\nTask {result['task_id']} FAILED: {result['error']}")
    
    # Summary
    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful
    print(f"\n{'='*60}")
    print(f"Completed: {successful}/{len(task_ids)} successful, {failed} failed")
    if successful > 0:
        total_time = sum(r['time'] for r in results if r['success'])
        print(f"Total compute time: {total_time:.2f}s")
