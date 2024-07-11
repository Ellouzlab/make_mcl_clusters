from scripts.utils import running_message, read_fasta, write_fasta
from scripts.analyze import analyze
from subprocess import run
from math import log10
import pandas as pd
import os

@running_message
def diamond(input, database, tsv_path, bitscore, threads, sensitivity=1):
    sensitivity_types = {
        1: "",
        2: "--more-sensitive",
        3: "--very-sensitive",
        4: "--ultra-sensitive"
    }
    if not sensitivity in sensitivity_types:
        raise ValueError("Sensitivity should be between 1 and 4")
    
    if not os.path.exists(tsv_path):
        cmd = f"diamond blastp {sensitivity_types[sensitivity]} -q {input} -d {database} --min-score {bitscore} -o {tsv_path} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --threads {threads}"
        run(cmd, shell=True)
    else: 
        print("TSV file already exists, using existing file")

@running_message
def makedb(input, database):
    if not os.path.exists(f"{database}.dmnd"):
        cmd=f"diamond makedb --in {input} -d {database}"
        run(cmd, shell=True)
    else:
        print("Database already exists, using existing database")

@running_message 
def mcl(input, output, inflation, threads):
    if not os.path.exists(output):
        cmd = f"mcl {input} --abc -I {inflation} -o {output} -te {threads}"
        run(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")

@running_message
def parsing_clusters(mcl_output):
    with open(mcl_output) as handle:
        cluster_lines = handle.readlines()
    
    id=1
    clust_dict = {}
    for line in cluster_lines:
        if "\t" in line:
            clust_dict[id] = line.strip().split("\t")
            id+=1
            
    gene_dict = {}
    for clust_num in clust_dict:
        for gene in clust_dict[clust_num]:
            gene_dict[gene] = clust_num
            
    
    return gene_dict


def cluster(input: str, outdir: str, threads: int)->None:
    '''
    Cluster proteins using MCL
    :param input: input fasta file
    :param outdir: output directory
    :param threads: number of threads
    :param perform_analysis: whether to perform the analysis step
    :return: None
    '''

    os.makedirs(outdir, exist_ok=True)
    
    # Make database for diamond
    database_path = f"{outdir}/database"
    makedb(input, database_path)
    
    # Run diamond
    tsv_path = f"{outdir}/diamond.tsv"
    diamond(input, database_path, tsv_path, 50, threads, 4)
    
    # Read pd dataframe from tsv
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
    df = pd.read_csv(tsv_path, sep="\t", names=columns)
    
    # Prepare for MCL
    df['evalue'] = df['evalue'].apply(lambda x: x if x > 0 else 1e-300)
    df["neglogeval"] = df['evalue'].apply(lambda x: -log10(x))
    
    mcl_input = f"{outdir}/mcl_input.tsv"
    mcl_df = df[["qseqid", "sseqid", "neglogeval"]]
    mcl_df.to_csv(mcl_input, sep='\t', header=False, index=False)
    
    # Cluster using MCL
    mcl(mcl_input, f"{outdir}/mcl_output", 1.3, threads)
    gene_dict = parsing_clusters(f"{outdir}/mcl_output")
    
    rep_gene_score_dict=analyze(df, gene_dict, outdir)
    
    # Write representative genes to a fasta file
    record_list = read_fasta(input)
    rep_gene_list = []
    for rep_gene in rep_gene_score_dict:
        for record in record_list:
            if record.id == rep_gene:
                rep_gene_list.append(record)
                break
    
    for record in rep_gene_list:
        record.description = f"{gene_dict[record.id]}"

    write_fasta(f"{outdir}/representative_genes.fasta", rep_gene_list)
