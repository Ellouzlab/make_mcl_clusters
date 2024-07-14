from scripts.cluster.cluster import cluster
from scripts.network.network import network
from scripts.arguments import arguments

def main():
    args = arguments()
    
    subcommands = ['cluster', 'network']
    
    if args.command == 'cluster':
        cluster(
            input=args.input, 
            outdir=args.outdir, 
            threads=args.threads,
            sensitivity=args.sensitivity)
    elif args.command == 'network':
        network(
            input=args.input, 
            map=args.map, 
            threads=args.threads, 
            outdir=args.outdir)
    else:
        print(f"No valid subcommand provided please specify args {subcommands}. Use --help for more information.")

    