import argparse
from multiprocessing import cpu_count

def arguments():
    version = "1.0"
    
    args = argparse.ArgumentParser(description="Make MCL cluster using Diamond")
    args.add_argument('-v', '--version', action='version', version=f'%(prog)s {version}')
    subparsers = args.add_subparsers(dest='command', help='sub-command help')
    
    cluster_parser = subparsers.add_parser('cluster', help='build cluster representative database')
    cluster_parser.add_argument("-i", "--input", help="Input protein fasta", required=True)
    cluster_parser.add_argument("-o", "--outdir", help="Output directory", default="output_cluster")
    cluster_parser.add_argument("-t", "--threads", help="Number of threads", default=cpu_count(), type=int)
    cluster_parser.add_argument("-s", "--sensitivity",
        help=(
            "Sensitivity of the clustering: "
            "0 - Faster, "
            "1 - Fast (default), "
            "2 - sensitive, "
            "3 - More sensitive, "
            "4 - Very sensitive, "
            "5 - Ultra sensitive"
        ),
        default=1,
        type=int,
        choices=[0, 1, 2, 3, 4, 5])
    cluster_parser.add_argument("--mmseqs", help="Use mmseqs instead of diamond", action="store_true")
    cluster_parser.add_argument("--mmseqs_sensitivity", help="Sensitivity of mmseqs clustering", default=7.5, type=float)
    
    network_parser = subparsers.add_parser('network', help='create a network from clusters')
    network_parser.add_argument("-i", "--input", help="Output directory from the 'cluster' module", required=True)
    network_parser.add_argument("-r", '--reference_db', help="TSV mapping file with first column 'assembly_name' being nodes, and second, 'protein_ids' being protein ids in the initial protein file", required=True)
    network_parser.add_argument("-t", "--threads", help="Number of threads", default=cpu_count(), type=int)
    network_parser.add_argument("-o", "--outdir", help="Output directory", default="output_network")
    network_parser.add_argument("--type", default='nucl', choices=['nucl', 'prot'], help="reference database type")
    
    gene_share_parser = subparsers.add_parser('gene_share', help='calculate gene sharing network')
    gene_share_parser.add_argument("-i", "--input", help="Output directory from the 'cluster' module", required=True)
    gene_share_parser.add_argument("-m", "--mapping", help="TSV mapping file with first column 'assembly_name' being nodes, and second, 'protein_ids' being protein ids in the initial protein file", default='None')
    gene_share_parser.add_argument("-t", "--threads", help="Number of threads", default=cpu_count(), type=int)
    gene_share_parser.add_argument("-o", "--outdir", help="Output directory", default="output_gene_share")
    gene_share_parser.add_argument('--gen_mapping_file', help = 'only use if your protein name is node_id_{some number} made for the developer', action='store_true', default=False)

    arguments = args.parse_args()

    #Throw errors for bad choices
    if arguments.command == "cluster":
        if arguments.mmseqs and arguments.sensitivity != 1:
            args.error("--mmseqs flag is incompatible with --sensitivity flag. Please use the default sensitivity (1) when using --mmseqs. To alter mmseqs sensitivity, use --mmseqs_sensitivity flag.")
            
        if not arguments.mmseqs and arguments.mmseqs_sensitivity != 7.5:
            args.error("--mmseqs_sensitivity flag is only compatible with --mmseqs flag.")
    
    if arguments.command == "gene_share":
        if arguments.gen_mapping_file and arguments.mapping != 'None':
            args.error("--gen_mapping_file flag is only compatible when mapping file not provided flag.")
            
        if not arguments.gen_mapping_file and arguments.mapping == 'None':
            args.error("Please provide a mapping file or use --gen_mapping_file flag to generate a mapping if node_id_{some number} = protein_id.")
    
    return arguments

if __name__ == "__main__":
    args = arguments()
