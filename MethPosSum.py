import os
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

def argparser():
    parser = ArgumentParser(
        description="Produce tsv file of methylated positions based on MethCalling output",
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=True
    )
    # Required arguments
    parser.add_argument("--input_file",
                        help='input file comes from methcalling results.',
                        required=True)
    parser.add_argument("--out_file",
                        help='filename to output summary results.',
                        required=True)

    # Optional arguments
    parser.add_argument("--min_read_depth",
                        help='minimum coverage of position to determine methylation.',
                        default=1, type=int)
    parser.add_argument('--mod_threshold',
                        help='minimum fraction of observations with probability of methylation >=50% at a position to include in report',
                        default=0.5, type=float)
    return parser


def check_thresh(locus_list, depth_thresh):

    if len(locus_list) >= depth_thresh:
        return True
    else:
        return False

def aggregate_by_pos(meth_file,output_file,min_read_depth,mod_threshold):

    pos_dict = {}
    for line in open(meth_file,'r'):
        if not line.startswith("contig"):
            try:
                contig, read, pos, strand, kmer, prob, label = tuple(line.strip("\n").split(','))
            except:
                print("number of columns in input file is wrong!")

            if kmer[int(len(kmer)/2)] != 'A':
                print("error kmer")

            if (contig,pos,strand,kmer) not in pos_dict:
                pos_dict[(contig,pos,strand,kmer)] = []

            if label == "1":
                pos_dict[(contig,pos,strand,kmer)].append(1)
            else:
                pos_dict[(contig,pos,strand,kmer)].append(0)

    # print(pos_dict)

    outfi = open(output_file, 'w')
    outfi.write("contig\tposition\tstrand\tref_11mer_sequence\tdepth_coverage\tfraction\tlabel\n")
    for locus in pos_dict.keys():
        a = check_thresh(pos_dict[locus], min_read_depth)
        if a:
            frac = np.mean(pos_dict[locus])
            if frac >= mod_threshold:
                meth = "m6A"
            else:
                meth = "A"
            coverage = len(pos_dict[locus])
            outfi.write('\t'.join(list(locus)[:])+'\t'+str(coverage)+'\t'+str(frac)+'\t'+meth+'\n')

    outfi.close()

def main(args):

    assert os.path.isfile(args.input_file), 'file not found at ' + args.input_file
    aggregate_by_pos(args.input_file, args.out_file, args.min_read_depth, args.mod_threshold)


if __name__=='__main__':
    # args = argparser().parse_args(["--input_file", "test_data/results.per_read_methylation", "--out_file", "test_data/results.per_sites_methylation.tsv"])
    args = argparser().parse_args()
    main(args)