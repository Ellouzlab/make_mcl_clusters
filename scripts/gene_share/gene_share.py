from scripts.utils import read_fasta, read_lines, pd_read_csv
import pandas as pd
import os
import numpy as np
from scipy.stats import hypergeom
from scipy.sparse import csr_matrix
from tqdm import tqdm

def calculate_gs_input(input_fasta, mcl_output, rep_fasta):
    print('Prepping diamond cluster output for gene share analysis')
    rep_records_dict = {record.id: record for record in read_fasta(rep_fasta)}

    cluster_lines = read_lines(mcl_output)
    all_genes = [record.id for record in read_fasta(input_fasta)]

    gene_dict = {}
    clu_num = 1
    for line in cluster_lines:
        for gene in line.split('\t'):
            gene_dict[gene.strip()] = clu_num
        clu_num += 1

    clust_dict = {}
    for gene in gene_dict:
        if gene in rep_records_dict:
            clust_dict[gene_dict[gene]] = gene

    for gene in all_genes:
        if gene not in gene_dict.keys():
            gene_dict[gene] = clu_num
            clu_num += 1

    for gene in gene_dict:
        if gene_dict[gene] not in clust_dict.keys():
            clust_dict[gene_dict[gene]] = gene

    df = pd.DataFrame({
        "Rep": [clust_dict[clust] for clust in gene_dict.values()],
        "Gene": list(gene_dict.keys())
    })
    return all_genes, df

def calculate_adjacency_matrix(presence_absence_matrix):
    phages = presence_absence_matrix.index
    total_proteins = presence_absence_matrix.shape[1]
    epsilon = 1e-100  # Small value to avoid log(0)
    
    # Convert presence-absence matrix to sparse matrix
    pamatrix = csr_matrix(presence_absence_matrix.values)

    print('Calculating shared matrix (dot product)')
    shared_matrix = pamatrix.dot(pamatrix.T).toarray()
    print('Done calculating shared matrix')

    # Calculate number of proteins for each phage
    a = pamatrix.sum(axis=1).A1  # Convert to 1D array
    b = pamatrix.sum(axis=1).A1  # Convert to 1D array

    # Vectorized calculation of hypergeometric p-values
    M = total_proteins
    N = a[:, None]  # Column vector
    K = b[None, :]  # Row vector

    pval_matrix = np.zeros_like(shared_matrix, dtype=float)
    for k in tqdm(range(shared_matrix.max() + 1)):
        pval_matrix += (shared_matrix >= k) * hypergeom.sf(k - 1, M, N, K)

    adjacency_matrix = -np.log10(pval_matrix + epsilon)
    adjacency_matrix_df = pd.DataFrame(adjacency_matrix, index=phages, columns=phages)

    return adjacency_matrix_df
def gene_share(input, mapping, threads, outdir, gen_map):
    if os.path.exists(f"{input}/cluster_output.tsv"):
        rep_df = pd_read_csv(f"{input}/cluster_output.tsv", sep='\t', names=['Rep', 'Gene'])
        all_genes = [record.id for record in read_fasta(f"{input}/input.fasta")]
    elif os.path.exists(f"{input}/mcl_output.txt") and os.path.exists(f"{input}/representative_genes.fasta"):
        all_genes, rep_df = calculate_gs_input(f"{input}/input.fasta", f"{input}/mcl_output.txt", f"{input}/representative_genes.fasta")
    else:
        print('No valid input files found, please run cluster command first')
        raise FileNotFoundError

    os.makedirs(outdir, exist_ok=True)

    if gen_map:
        print("Generating presence absence matrix")
        
        gene_to_organism = {gene: '_'.join(gene.split('_')[:-1]) for gene in all_genes}
        
        rep_df['Node'] = rep_df['Gene'].apply(lambda x: gene_to_organism[x])
        presence_absence_matrix = (rep_df.groupby(['Node', 'Rep'])
                                           .size()
                                           .unstack(fill_value=0))
        
        
        print("Calculating adjacency matrix")
        
        # Calculate adjacency matrix with vectorized operations
        adjacency_matrix = calculate_adjacency_matrix(presence_absence_matrix)
        adjacency_matrix_file = f"{outdir}/adjacency_matrix.csv"
        adjacency_matrix.to_csv(adjacency_matrix_file)
        print(f"Adjacency matrix saved to {adjacency_matrix_file}")
    else:
        print('Gene mapping not provided. Skipping the creation of the presence-absence matrix.')

