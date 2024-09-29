# README

This program is designed to call 6mA DNA methylation from nanopore data using the deeplearning (For prokaryotes).

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
* multi_to_single_fast5
* nanopolish (https://github.com/jts/nanopolish) or f5c (recommend)
* minimap2
* samtools
* slow5tools (option)
* nanopore sequencing data (fastq format + fast5 to run nanopolish, basecalled using Guppy or another basecaller that saves event data)
* a reference sequence file (fasta)

# Installation using conda
    git clone https://github.com/caojiabao/methCaller.git
    conda env create -f environment.yml
    conda activate methCaller
    
then add softlinks for MethCalling.py, MethPosSum.py and PrepareData.py to path or run from MethCalling directory

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

# Pipeline for 6mA methylation detection from R9.4.1 data

#### 1. extract template strand reads from fast5 files.
~~~
multi_to_single_fast5 -i <multiple_fast5_dir> -s <single_fast5_dir> -t 20 --recursive
f5c_x86_64_linux index -d <single_fast5_dir> <reads>.fastq -t 8 -s sequencing_summary.txt
~~~
#### 2. align fastq reads to reference
~~~
minimap2 -ax map-ont <reference>.fasta <reads>.fastq -t 8 > <reads>.sam
samtools sort -T tmp -@ 8 -o <reads>_sorted.bam <reads>.sam
samtools index <reads>_sorted.bam
~~~
#### 3. event align using f5c (recommend) or nanopolish
~~~
f5c_x86_64_linux eventalign -r <reads>.fastq -b <reads>_sorted.bam -g <reference>.fasta -t 8 --iop 8 --scale-events --signal-index > <reads>_eventalign.tsv
~~~
#### 4. 6mA methylation prediction using methCaller
~~~
python PrepareData.py --eventalign <reads>_eventalign.tsv --fasta_file <reference>.fasta --out_dir <results_dir> --n_processes 8 --chunk_size 40960

python MethCalling.py --input_file <results_dir>/data.info --out_dir <results_dir> --model_dir models/biLSTM --model bilstm --n_processes 8 --chunk_size 4096

python MethPosSum.py --input_file <results_dir>/results.per_read_methylation --out_file <results_dir>/results.per_sites_methylation.tsv
~~~

# Explanation of results
At the end of the program run, five files are generated, of which `results.per_sites_methylation.tsv` is the final result of whether or not methylation occurs at adenine sites in the reference genome.
##### The file contains a total of 7 columns:
1. contig column: the name of the chromosome in the reference genome;
2. position column: position of adenine in the chromosome (0-based); the
3. strand column: positive and negative strands;
4. ref_11mer_sequence column: indicates the 11-mer sequence at the position;
5. depth_coverage column: the sequencing depth of the position;
6. fraction column; proportion of methylated reads, fraction >= 0.5 is considered methylated;
7. label column; methylated or not

# Test data
The reference genome sequence and eventalign results have been saved in the test_data directory.
##### 1. extraction feature from eventalign results
~~~ 
python PrepareData.py --eventalign test_data/test_eventalign.tsv --fasta_file test_data/test_contigs.fa --out_dir test_data --n_processes 8 --chunk_size 256
~~~
##### 2. prediction of 6mA methylation using methCaller
~~~
python MethCalling.py --input_file test_data/data.info --out_dir test_data --model_dir models/biLSTM --model bilstm --n_processes 8 --chunk_size 256
~~~
##### 3. summarize 6mA methylation prediction results
~~~
python MethPosSum.py --input_file test_data/results.per_read_methylation --out_file test_data/results.per_sites_methylation.tsv
~~~
