import os
import torch
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from Utils.model import BiLSTM_Attention, device, load_checkpoint
from Utils.predict import run_predict

def argparser():
    parser = ArgumentParser(
        description="methylation prediction of per read mapped to reference genome using pretrained models",
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=True
    )
    # Required arguments
    parser.add_argument("--input_file",
                        help='data.info file comes from forward step PrepareData.',
                        required=True)
    parser.add_argument("--out_dir",
                        help='directory to output prediction results.',
                        required=True)
    parser.add_argument("--model_dir",
                        help="directory to pretrained models",
                        required=True)

    # Optional arguments
    parser.add_argument("--model",
                        help="choice of deep learning model",
                        default="bilstm", type=str)
    parser.add_argument("--n_processes",
                        help='number of processes to run.',
                        default=1, type=int)
    parser.add_argument('--chunk_size',
                        help='number of lines from data.info for processing.',
                        default=256, type=int)
    parser.add_argument("--seed",
                        help='random seed for sampling.',
                        default=0, type=int)
    parser.add_argument("--read_proba_threshold",
                        help='default probability threshold for a read to be considered modified',
                        default=0.5, type=float)
    return parser


def main(args):

    input_file = args.input_file
    out_dir = args.out_dir

    torch.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)
    np.random.seed(args.seed)

    #biLSTM model
    if args.model == "bilstm":
        model = BiLSTM_Attention().to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        model_path = args.model_dir + "/bilstm_attention_1X7_917.pt"
        load_checkpoint(model_path, model, optimizer)
    else:
        print("model can not empty")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # with open(os.path.join(args.out_dir, "data.site_proba.csv"),'w', encoding='utf-8') as f:
    #     f.write('contig_id,contig_position,n_reads,probability_modified,kmer,mod_ratio\n')
    # with open(os.path.join(args.out_dir, "data.read_proba.csv"), 'w', encoding='utf-8') as g:
    #     g.write('contig_id,contig_position,read_index,probability_modified\n')
    run_predict(input_file, out_dir, args, model)

if __name__=='__main__':
    # args = argparser().parse_args(["--input_file", "test_data/data.info", "--out_dir", "test_data", "--model_dir", "models/biLSTM","--model", "bilstm", "--n_processes", "2"])
    torch.multiprocessing.set_start_method("spawn")
    args = argparser().parse_args()
    main(args)


