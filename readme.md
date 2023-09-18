# README

This program is designed to call 6mA DNA methylation from nanopore data using the deeplearning.

# Dependencies/requirements

**python packages**
* python3.9
* torch
* numpy
* biopython
* pandas
* re
* os
* multiprocessing

**other**
* nanopolish (https://github.com/jts/nanopolish) or f5c (recommend)
* minimap2
* samtools
* slow5tools
* nanopore sequencing data (fastq format + fast5 to run nanopolish, basecalled using Guppy or another basecaller that saves event data)
* a reference sequence file (fasta)

# Installation
Add softlinks for MethCalling.py, MethPosSum.py and PrepareData.py to path or run from MethCalling directory

# Options
    usage: MethCalling.py [-h] --input_file INPUT_FILE 
            --out_dir OUT_DIR [--n_processes N_PROCESSES] 
            [--chunk_size CHUNK_SIZE] [--seed SEED] 
            [--read_proba_threshold READ_PROBA_THRESHOLD]


**optional arguments:**
  
    --h, --help            show this help message and exit
    --input_file INPUT_FILE
                        data.info file comes from forward step PrepareData. (default: None)
    --out_dir OUT_DIR     directory to output prediction results. (default: None)
    --n_processes N_PROCESSES
                        number of processes to run. (default: 1)
    --chunk_size CHUNK_SIZE
                        number of lines from data.info for processing. (default: 256)
    --seed SEED           random seed for sampling. (default: 0)
    --read_proba_threshold READ_PROBA_THRESHOLD
                        default probability threshold for a read to be considered modified (default: 0.5)

