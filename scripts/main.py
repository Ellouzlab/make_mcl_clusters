from scripts.cluster.cluster import cluster
from scripts.cluster.mmseqs_cluster import mmseqs_cluster
from scripts.network.network import network
from scripts.gene_share.gene_share import gene_share
from scripts.arguments import arguments


def main():
    args = arguments()
    subcommands = ['cluster', 'network', 'gene_share']
    

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
            type=args.type)
    elif args.command == 'gene_share':
        gene_share(
            input=args.input, 
            mapping=args.mapping, 
            threads=args.threads, 
            outdir=args.outdir,
            gen_map = args.gen_mapping_file
        )
    else:
        print(f"No valid subcommand provided please specify: {subcommands}. Use --help for more information.")

    