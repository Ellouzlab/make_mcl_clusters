import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scripts.network.mmseqs_network import mmseqs_search
from scripts.cluster.mmseqs_cluster import mmseqs_makedb
from scripts.utils import pd_read_csv, read_fasta

def filter_hits(df, bitscore_threshold):
    return df[(df['bitscore'] >= bitscore_threshold) & (~df['query'].str.contains('Pantoea'))]


def calculate_proximity_counts(viral_host_search_res, protein_distance):
    alone_counts = defaultdict(int)
    total_counts = defaultdict(int)
    
    for query, group in viral_host_search_res.groupby('query'):
        group = group.sort_values('tstart')
        prev_end = None
        for _, row in group.iterrows():
            total_counts[query] += 1
            if prev_end is None or row['tstart'] - prev_end > protein_distance:
                alone_counts[query] += 1
            prev_end = row['tend']
    
    return alone_counts, total_counts

def categorize_viral_pcs(viral_host_search_res, viral_go_search_res, bitscore_threshold):
    with_go_metabolic_hits = set(viral_go_search_res[viral_go_search_res['bitscore'] >= bitscore_threshold]['query'])
    without_go_metabolic_hits = set(viral_host_search_res['query']) - with_go_metabolic_hits
    
    return with_go_metabolic_hits, without_go_metabolic_hits

def plot_histogram(alone_counts, total_counts, with_go_metabolic_hits, without_go_metabolic_hits):
    # Prepare data for histogram
    ratios_with_go = {pc: alone_counts[pc] / total_counts[pc] for pc in with_go_metabolic_hits if total_counts[pc] > 10}
    ratios_without_go = {pc: alone_counts[pc] / total_counts[pc] for pc in without_go_metabolic_hits if total_counts[pc] > 10}
    
    # Data for plotting
    data_with_go = list(ratios_with_go.values())
    data_without_go = list(ratios_without_go.values())

    plt.figure(figsize=(10, 6))
    plt.hist(data_with_go, bins=20, alpha=0.75, color='blue', label='With GO Metabolic Hits')
    plt.hist(data_without_go, bins=20, alpha=0.75, color='red', label='Without GO Metabolic Hits')
    plt.xlabel('Alone / Total Hits Ratio')
    plt.ylabel('Frequency')
    plt.title('Histogram of Alone/Total Hits Ratios for Viral PCs')
    plt.legend()
    plt.grid(True)
    plt.show()

def process_mmseqs(search_outdir):
    columns = ["query", "target", "pident", "alnlen", "mismatch", "numgapopen", "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]
    df = pd_read_csv(f"{search_outdir}/network.m8", sep="\t", names=columns)
    return df

def amg(viral_pcs, go_metabolic, host_genomes, protein_distance, threads, outdir):
    os.makedirs(outdir, exist_ok=True)
    bitscore_threshold = 10
    protein_distance = 100000

    db = f"{outdir}/db"
    search_results = f"{outdir}/search_results"

    viral_pcs_db = f"{db}/viral_pcs_db"
    go_metabolic_db = f"{db}/go_metabolic_db"
    host_genomes_db = f"{db}/host_genomes_db"
    os.makedirs(viral_pcs_db, exist_ok=True)
    os.makedirs(go_metabolic_db, exist_ok=True)
    os.makedirs(host_genomes_db, exist_ok=True)

    mmseqs_makedb(viral_pcs, viral_pcs_db)
    mmseqs_makedb(go_metabolic, go_metabolic_db)
    mmseqs_makedb(host_genomes, host_genomes_db)

    viral_metabolic_outdir = f"{search_results}/viral_metabolic"
    if not os.path.exists(viral_metabolic_outdir):
        mmseqs_search(viral_pcs_db, go_metabolic_db, viral_metabolic_outdir)
    
    host_viral_outdir = f"{search_results}/host_viral"
    if not os.path.exists(host_viral_outdir):
        mmseqs_search(viral_pcs_db, host_genomes_db, host_viral_outdir, addition="--search-type 4")
    
    viral_go_search_res = process_mmseqs(viral_metabolic_outdir)
    viral_host_search_res = process_mmseqs(host_viral_outdir)
    viral_go_search_res = filter_hits(viral_go_search_res, bitscore_threshold)
    viral_host_search_res = filter_hits(viral_host_search_res, bitscore_threshold)

    alone_counts, total_counts = calculate_proximity_counts(viral_host_search_res, protein_distance)
    with_go_metabolic_hits, without_go_metabolic_hits = categorize_viral_pcs(viral_host_search_res, viral_go_search_res, bitscore_threshold)

    plot_histogram(alone_counts, total_counts, with_go_metabolic_hits, without_go_metabolic_hits)

# Configuration and call for the function here if necessary
