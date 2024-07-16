from scripts.utils import read_fasta, read_lines, pd_read_csv
import pandas as pd
import os
from tqdm import tqdm

def calculate_gs_input(input_fasta ,mcl_output, rep_fasta):
    print('prepping diamond cluster output for gene share analysis')
    rep_records_dict = {record.id :record for record in read_fasta(rep_fasta)}

    cluster_lines = read_lines(mcl_output)
    all_genes = [record.id for record in read_fasta(input_fasta)]

    gene_dict = {}
    clu_num=1
    for line in cluster_lines:
        for gene in line.split('\t'):
            gene_dict[gene.strip()] = clu_num
        clu_num+=1
    
    clust_dict={}
    for gene in gene_dict:
        if gene in rep_records_dict:
            clust_dict[gene_dict[gene]] = gene

    for gene in all_genes:
        if not gene in gene_dict.keys():
            gene_dict[gene] = clu_num
            clu_num+=1

    for gene in gene_dict:
        if not gene_dict[gene] in clust_dict.keys():
            clust_dict[gene_dict[gene]] = gene
    
    df = pd.DataFrame({
        "Rep": [clust_dict[clust] for clust in gene_dict.values()],
        "Gene": list(gene_dict.keys())
    })
    return all_genes, df

def gene_share(input, mapping, threads, outdir, gen_map):
    
    if os.path.exists(f"{input}/cluster_output.tsv"):
        rep_df = pd_read_csv(f"{input}/cluster_output.tsv", sep='\t', names=['Rep', 'Gene'])
        all_genes = [record.id for record in read_fasta(f"{input}/input.fasta")]
    
    elif os.path.exists(f"{input}/mcl_output.txt") and os.path.exists(f"{input}/representative_genes.fasta"):
        all_genes, rep_df = calculate_gs_input(f"{input}/input.fasta",f"{input}/mcl_output.txt", f"{input}/representative_genes.fasta")
    
    else:
        print('No valid input files found, please run cluster command first')
        raise FileNotFoundError
    
    # Create mapping from genes to organisms
    if gen_map:
        gene_to_organism = {gene: '_'.join(gene.split('_')[:-1]) for gene in all_genes}
        
        # Create a DataFrame to hold the presence-absence matrix
        organisms = rep_df['Gene'].apply(lambda x: gene_to_organism[x]).unique()
        clusters = rep_df['Rep'].unique()
        presence_absence_matrix = pd.DataFrame(0, index=organisms, columns=clusters)
        
        # Fill the matrix with counts
        for _, row in tqdm(rep_df.iterrows(), total=rep_df.shape[0], desc="Processing genes"):
            organism = gene_to_organism[row['Gene']]
            cluster = row['Rep']
            presence_absence_matrix.loc[organism, cluster] += 1

        # Save the presence-absence matrix
        presence_absence_matrix.to_csv(os.path.join(outdir, 'presence_absence_matrix.csv'))
        print(f'Presence-absence matrix saved to {os.path.join(outdir, "presence_absence_matrix.csv")}')
    else:
        print('Gene mapping not provided. Skipping the creation of the presence-absence matrix.')
