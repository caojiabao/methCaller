import warnings
import os
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
import pandas as pd
from Utils.Utils import eventalign_index_parallel
from Utils.dataprepare import parallel_preprocess_tx

def argparser():
    parser = ArgumentParser(
        description="extraction of 11-mer sequence features from nanopolish eventalign results",
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=True
    )

    # Required arguments
    parser.add_argument('--eventalign',
                        help='eventalign filepath, the output from nanopolish.',
                        required=True)

    parser.add_argument('--fasta_file',
                        help='reference genomes.',
                        required=True)

    parser.add_argument('--out_dir',
                        help='output directory.',
                        required=True)


    # Optional arguments

    parser.add_argument('--n_processes',
                        help='number of processes to run.',
                        default=1, type=int)
    parser.add_argument('--chunk_size',
                        help='number of lines from nanopolish eventalign.txt for processing.',
                        default=10000, type=int)
    parser.add_argument('--readcount_min',
                        help='minimum read counts per gene',
                        default=1, type=int)

    return parser


def main(args):

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

    eventalign_index_parallel(args.eventalign, args.chunk_size, args.out_dir, args.n_processes)

    # For each read, combine multiple events aligned to the same positions,
    # the results from nanopolish eventalign, into a single event per position.
    parallel_preprocess_tx(args.eventalign, args.fasta_file, args.out_dir, args.n_processes, args.readcount_min)

if __name__ == '__main__':
    # args = argparser().parse_args(["--eventalign","test_data/test_eventalign.tsv","--fasta_file","test_data/test_contigs.fa","--out_dir","test_data"])
    args = argparser().parse_args()
    main(args)

