from src import __version__, __author__
from os import cpu_count
import argparse


def arguments():
    
    args=argparse.ArgumentParser(description="Make MCL cluster using Diamond")
    args.add_argument("-i", "--input", help="Input protein fasta", required=True)
    args.add_argument("-o", "--outdir", help="Output directory", default="output_mcl_cluster")
    args.add_argument("-t", "--threads", help="Number of threads", default=cpu_count(), type=int)
    args.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__} by {__author__}')
    return args.parse_args()
