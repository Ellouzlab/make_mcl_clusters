from scripts.network.dna_translator import translate
from scripts.network.diamond_network import diamond_network
from scripts.network.mmseqs_network import mmseqs_network
import os

def network(input, reference, threads, outdir, type, mmseqs):
    rep_seq_path = f"{input}/representative_genes.fasta"
    
    os.makedirs(outdir, exist_ok=True)
    
    if type == 'nucl':
        translated_path = f"{outdir}/translated_reference_db.fasta"
        translate(reference, translated_path, threads)
    else:
        translated_path = reference
    
    if mmseqs:
        mmseqs_network(
            query = rep_seq_path,
            outdir = outdir,
            threads = threads,
            reference_db = translated_path,
        )
    else:
        diamond_network(
            query = rep_seq_path,
            outdir = outdir,
            threads = threads,
            reference_db = translated_path,
        )
    
    
    