from scripts.utils import running_message, pd_read_csv, read_fasta, write_fasta
from subprocess import run
from tqdm import tqdm
import os
from scripts.mmseqs_utils import mmseqs_makedb, mmseqs_cluster_cmd, mmseqs_createtsv

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

@running_message
def mmseqs_cluster(input: str, outdir: str, threads: int, sensitivity: int)->None:
    os.makedirs(outdir, exist_ok=True)
    
    fasta_path = f"{outdir}/input.fasta"
    if not os.path.exists(fasta_path):
        run(f"cp {input} {fasta_path}", shell=True)
    
    db = mmseqs_makedb(input, outdir)
    cluster_output = mmseqs_cluster_cmd(db, outdir, threads, sensitivity)
    tsv_path = mmseqs_createtsv(db, cluster_output, outdir)
    
    processing_cluster(tsv_path, input, outdir)
    

    
    
        