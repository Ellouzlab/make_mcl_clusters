from scripts.utils import running_message, pd_read_csv, read_fasta, write_fasta
from subprocess import run
from tqdm import tqdm
import pandas as pd
import os

@running_message
def mmseqs_makedb(input, outdir):
    db = f"{outdir}/mmseqs_db"
    if not os.path.exists(f"{outdir}/mmseqs_db"):
        cmd = f"mmseqs createdb {input} {outdir}/mmseqs_db"
        run(cmd, shell=True)
    else:
        print("Database already exists, using existing database")
    return db

@running_message
def mmseqs(db, outdir, threads, sensitivity, coverage=0.5, min_id = 0.3 ,e_value=1e5):
    tmp = f"{outdir}/tmp"
    os.makedirs(tmp, exist_ok=True)
    
    cluster_output = f"{outdir}/cluster_output"
    if not os.path.exists(f"{cluster_output}.index"):
        cmd = f"mmseqs cluster -s {sensitivity} -c {coverage} --min-seq-id {min_id} -e {e_value} --threads {threads} {db} {cluster_output} {tmp}"
        run(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return cluster_output

@running_message
def mmseqs_createtsv(db, cluster_output, outdir):
    tsv_path = f"{outdir}/cluster_output.tsv"
    if not os.path.exists(tsv_path):
        cmd = f"mmseqs createtsv {db} {db} {cluster_output} {tsv_path}"
        run(cmd, shell=True)
    else:
        print("Output file already exists, using existing file")
    return tsv_path

@running_message
def processing_cluster(tsv_path, input, outdir):
    columns = ["rep", "gene"]
    df = pd_read_csv(tsv_path, sep="\t", names=columns)

    rep_gene_counts = df.groupby('rep')['gene'].nunique()

    reps_to_keep = rep_gene_counts.index.to_list()
    
    all_records = read_fasta(input)
    record_dict = {record.id: record for record in tqdm(all_records, desc="Creating record dictionary", unit=" Records")}
    
    rep_records=[]
    for rep in tqdm(reps_to_keep, desc="Matching records to fasta", unit=" Reps"):
        rep_records.append(record_dict[rep])
    
    write_fasta(f"{outdir}/representative_genes.fasta", rep_records)

def mmseqs_cluster(input: str, outdir: str, threads: int, sensitivity: int)->None:
    os.makedirs(outdir, exist_ok=True)
    
    fasta_path = f"{outdir}/input.fasta"
    if not os.path.exists(fasta_path):
        run(f"cp {input} {fasta_path}", shell=True)
    
    db = mmseqs_makedb(input, outdir)
    cluster_output = mmseqs(db, outdir, threads, sensitivity)
    tsv_path = mmseqs_createtsv(db, cluster_output, outdir)
    
    processing_cluster(tsv_path, input, outdir)
    

    
    
        