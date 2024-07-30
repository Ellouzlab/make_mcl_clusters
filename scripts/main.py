from scripts.cluster.cluster import cluster
from scripts.cluster.mmseqs_cluster import mmseqs_cluster
from scripts.network.network import network
from scripts.gene_share.gene_share import gene_share
from scripts.arguments import arguments
from scripts.amg.amg import amg
from scripts.utils import init_logging
from sys import exit
import os


def main():
    args = arguments()

    subcommands = ['cluster', 'network', 'gene_share', 'amg']
    
    if not args.command in subcommands:
        print(f"No valid subcommand provided please specify: {subcommands}. Use --help for more information.")
        exit(1)
    
    os.makedirs(args.outdir, exist_ok=True)
    log_filename = f"{args.outdir}/output.log"
    init_logging(log_filename)
    

    if args.command == 'cluster':
        if args.mmseqs:
            
            mmseqs_cluster(
                input=args.input, 
                outdir=args.outdir, 
                threads=args.threads, 
                sensitivity=args.mmseqs_sensitivity)
        else:
            cluster(
                input=args.input, 
                outdir=args.outdir, 
                threads=args.threads,
                sensitivity=args.sensitivity)
    elif args.command == 'network':
        network(
            input=args.input, 
            reference=args.reference_db, 
            threads=args.threads, 
            outdir=args.outdir,
            type=args.type,
            mmseqs=args.mmseqs)
    elif args.command == 'gene_share':
        gene_share(
            input=args.input, 
            mapping=args.mapping, 
            threads=args.threads, 
            outdir=args.outdir,
            gen_map = args.gen_mapping_file
        )
    elif args.command == 'amg':
        amg(
            viral_pcs=args.viral_pcs, 
            go_metabolic=args.go_metabolic, 
            host_genomes=args.host_genomes, 
            protein_distance=args.protein_distance, 
            threads=args.threads, 
            outdir=args.outdir
        )