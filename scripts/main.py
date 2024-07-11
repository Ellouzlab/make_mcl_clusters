from scripts.cluster import cluster
from scripts.arguments import arguments

def main():
    args = arguments()
    
    cluster(
        input=args.input, 
        outdir=args.outdir, 
        threads=args.threads)

    