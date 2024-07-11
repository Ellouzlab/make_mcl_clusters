from scripts.utils import running_message
from os import makedirs
from os.path import exists
from matplotlib import pyplot as plt
from tqdm import tqdm
import pandas as pd

def graphics(same_df, opposite_df, outdir):
    plt.figure(figsize=(10, 6))
    plt.hist(same_df['bitscore'], bins=5000, alpha=0.5, label='Same Cluster', color='blue')
    plt.hist(opposite_df['bitscore'], bins=5000, alpha=0.5, label='Opposite Cluster', color='red')
    plt.ylabel('Frequency')
    plt.title('Histogram of Bitscores')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(f"{outdir}/bitscore_histogram.png")
    plt.close()

def find_representative_gene(same_df, gene_dict):
    gene_points_df = same_df.groupby('qseqid')['bitscore'].sum().reset_index()
    gene_points_df['cluster'] = gene_points_df['qseqid'].map(gene_dict)
    gene_points_df = gene_points_df.dropna(subset=['cluster'])
    idx = gene_points_df.groupby('cluster')['bitscore'].idxmax()
    representative_genes = gene_points_df.loc[idx]
    representative_genes_dict = representative_genes.set_index('qseqid')['bitscore'].to_dict()

    return representative_genes_dict

@running_message
def analyze(df, gene_dict, outdir):
    def filter(row):
        try:
            if gene_dict[row["qseqid"]] == gene_dict[row["sseqid"]]:
                return True
            else:
                return False
        except:
            return False
    
    path_to_same_cluster = f"{outdir}/same_cluster.tsv"
    path_to_opposite_cluster = f"{outdir}/opposite_cluster.tsv"

    if exists(path_to_same_cluster) and exists(path_to_opposite_cluster):
        print("Output file already exists, using existing file")
        same_cluster = pd.read_csv(path_to_same_cluster, sep="\t")
        opposite_cluster = pd.read_csv(path_to_opposite_cluster, sep="\t")
    else:
        mask = df.apply(filter, axis=1)

        same_cluster = df[mask]
        opposite_cluster = df[~mask]

        same_cluster.to_csv(path_to_same_cluster, sep="\t", index=False)
        opposite_cluster.to_csv(path_to_opposite_cluster, sep="\t", index=False)
    
    new_dir = f"{outdir}/graphics"

    makedirs(new_dir, exist_ok=True)
    graphics(same_cluster, opposite_cluster, new_dir)
    
    return find_representative_gene(same_cluster, gene_dict)

    



    
