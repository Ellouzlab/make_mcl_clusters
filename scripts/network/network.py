from scripts.cluster.cluster import makedb, diamond
from scripts.network.dna_translator import translate
import os

def network(input, reference, threads, outdir, type):
    rep_seq_path = f"{input}/representative_genes.fasta"
    
    os.makedirs(outdir, exist_ok=True)
    db_path = f"{outdir}/reference_database"
    tsv_path = f"{outdir}/diamond_output.tsv"
    
    if type == 'nucl':
        translated_path = f"{outdir}/translated_reference_db.fasta"
        translate(reference, translated_path, threads)
        makedb(translated_path, db_path)
    else:
        makedb(reference, db_path)
    
    diamond(input = rep_seq_path,
            database = db_path,
            tsv_path = tsv_path,
            bitscore = 50,
            threads = threads,
            sensitivity = 5)