from src import __author__, __version__
from src.cluster import cluster
from src.arguments import arguments

def main():
    args = arguments()
    
    cluster(
        input=args.input, 
        outdir=args.outdir, 
        threads=args.threads)

    