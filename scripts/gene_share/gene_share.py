from scripts.utils import read_fasta, read_lines, pd_read_csv
import pandas as pd
import os

def calculate_gs_input(mcl_output, rep_fasta):
    print('prepping diamond cluster output for gene share analysis')
    rep_records_dict = {record.id :record for record in read_fasta(rep_fasta)}

    cluster_lines = read_lines(mcl_output)

    all_genes = []
    for line in cluster_lines:
        for gene in line.split('\t'):
            all_genes.append(gene)
    print(len(all_genes))

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

    for gene in gene_dict:
        if not gene_dict[gene] in clust_dict.keys():
            clust_dict[gene_dict[gene]] = gene
    
    df = pd.DataFrame({
        "Rep": [clust_dict[clust] for clust in gene_dict.values()],
        "Gene": list(gene_dict.keys())
    })
    return df



def gene_share(input, mapping, threads, outdir, gen_map):
    
    if os.path.exists(f"{input}/cluster_output.tsv"):
        rep_df = pd_read_csv(f"{input}/cluster_output.tsv", sep='\t', names=['Rep', 'Gene'])
    
    elif os.path.exists(f"{input}/mcl_output.txt") and os.path.exists(f"{input}/representative_genes.fasta"):
        rep_df = calculate_gs_input(f"{input}/mcl_output.txt", f"{input}/representative_genes.fasta")
    
    else:
        print('No valid input files found, please run cluster command first')
        raise FileNotFoundError
    
    print(rep_df)