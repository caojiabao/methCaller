import os,re
import multiprocessing
import numpy as np
import pandas as pd
from Bio.SeqIO import parse
from .helper import Consumer,end_queue

def eventalign_index(eventalign_result: pd.DataFrame, pos_start: int, out_paths, locks):
    r"""
    Function to index the position of a specific read features within eventalign.tsv
    :param eventalign_result: A pd.DataFrame object containing a portion of eventalign.tsv to be indexed
    :param pos_start: An index position within eventalign.tsv that corresponds to the start of the eventalign_result portion within eventalign.tsv file
    :param out_paths: A dictionary containing filepath for all the output files produced by the index function
    :param locks: A lock object from multiprocessing library that ensures only one process write to the output file at any given time
    :return: none
    """

    eventalign_result = eventalign_result.set_index(['contig','read_index'])
    pos_end=pos_start
    with locks['index'], open(out_paths['index'],'a', encoding='utf-8') as f_index:
        for _index in list(dict.fromkeys(eventalign_result.index)):
            contig,read_index = _index
            pos_end += eventalign_result.loc[_index]['line_length'].sum()
            f_index.write('%s,%d,%d,%d\n' %(contig,read_index,pos_start,pos_end))
            pos_start = pos_end

def eventalign_index_parallel(eventalign_filepath: str, chunk_size:int , out_dir:str, n_processes:int):
    r"""
    Function to index every read within eventalign.tsv file for faster access later
    :param eventalign_filepath:  String filepath to the eventalign.tsv file
    :param chunk_size: Chunksize argument for pd.read_csv function
    :param out_dir: String filepath to the output directory of the indexing function
    :param n_processes: Number of processes used for indexing
    :return: none
    """
    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['index']:
        out_paths[out_filetype] = os.path.join(out_dir,'eventalign.%s'  %out_filetype)
        locks[out_filetype] = multiprocessing.Lock()

    with open(out_paths['index'],'w', encoding='utf-8') as f:
        f.write('contig,read_index,pos_start,pos_end\n') # header

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.p
    consumers = [Consumer(task_queue=task_queue,task_function=eventalign_index,locks=locks) for i in range(n_processes)]
    for process in consumers:
        process.start()

    ## Load tasks into task_queue. A task is eventalign information of one read.
    eventalign_file = open(eventalign_filepath,'r', encoding='utf-8')
    pos_start = len(eventalign_file.readline()) #remove header
    chunk_split = None
    index_features = ['contig','read_index','line_length']
    for chunk in pd.read_csv(eventalign_filepath, chunksize=chunk_size,sep='\t'):
        chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
        chunk_concat = pd.concat([chunk_split,chunk_complete])
        chunk_concat_size = len(chunk_concat.index)
        ## read the file at where it left off because the file is opened once ##
        lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
        chunk_concat.loc[:, 'line_length'] = np.array(lines)
        task_queue.put((chunk_concat[index_features], pos_start, out_paths))
        pos_start += sum(lines)
        chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']].copy()

    ## the loop above leaves off w/o adding the last read_index to eventalign.index
    chunk_split_size = len(chunk_split.index)
    lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
    chunk_split.loc[:, 'line_length'] = np.array(lines)
    task_queue.put((chunk_split[index_features], pos_start, out_paths))

    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue,n_processes)

    # Wait for all the tasks to finish.
    task_queue.join()


def position_filtering(pos_list:list, sequence:str):
    r"""
    Function to filter positions located in start(n=5) and end(n=5) of reference
    """
    filtered_li = list()
    for i in pos_list:
        pos = i.span()[0]
        if pos >= 5 and pos <= (len(sequence) - 5):
                filtered_li.append(pos)
    return filtered_li

def fasta_to_position(fasta:str):
    r"""
    Function to search position of A base in reference genomes
    :param fasta: reference genomes
    :return: dict, contains locations of A base in reference genomes
    """
    chrome_positions = {}
    file = open(fasta)
    for record in parse(file, "fasta"):
        contig_seq = str(record.seq)
        contig_id = record.id
        li_A = re.finditer("A",contig_seq)
        li_A_filter = position_filtering(li_A,contig_seq)
        li_T = re.finditer("T", contig_seq)
        li_T_filter = position_filtering(li_T, contig_seq)
        base_positions = {}
        base_positions["A"] = li_A_filter
        base_positions["T"] = li_T_filter
        chrome_positions[contig_id] = base_positions
    return chrome_positions #dict:{A:[1,3],T:[1,5]}

def partions_continuous(pos_arr: np.ndarray, positions:list):
    r"""
    Function to extract indices of mapped 11-mers segmentions (filtering eventalign pass)
    :param pos_arr: reference position of read mapping from eventalin.tsv
    :param positions: location of A base in reference genomes
    :return: list
    """
    seg_indices = list()
    pos_arr = pos_arr.tolist()
    intersection = list(set(pos_arr) & set(positions))
    intersection.sort()
    for pos in intersection:
        indices = pos_arr.index(pos)
        start = indices - 5
        end = indices
        if start >= 0 and end <= (len(pos_arr) - 1):
            if pos_arr[end]-pos_arr[start]==5:
                t = (start,end+1)
                seg_indices.append(t)

    return seg_indices

def concensus_kmer(kmer_arr:np.ndarray,seg_indices:list):
    ref_11mers = []
    kmer_arr = kmer_arr.tolist()
    for i in seg_indices:
        kmer_list = kmer_arr[i[0]:i[1]]
        kmer_ref = kmer_list[0][0]
        for seq in kmer_list[1:]:
            kmer_ref += seq[0][-1]
        ref_11mers.append(kmer_ref)
    return ref_11mers

def concensus_model_kmer(model_kmer_arr:np.ndarray,seg_indices:list):
    ref_11mers = []
    kmer_arr = model_kmer_arr.tolist()
    for i in seg_indices:
        kmer_list = kmer_arr[i[0]:i[1]]
        kmer_ref = kmer_list[-1][0][:-1]
        for seq in kmer_list[::-1]:
            kmer_ref += seq[0][-1]
        ref_11mers.append(kmer_ref)
    return ref_11mers

def convert_event(event_array:np.recarray,index:list):
    arr_li = []
    for i in index:
        start = i[0]-1
        end = i[1]-1
        for x in range(end,start,-1):
            events = event_array[x]
            arr_li.append(events)
    return arr_li

def filter_events(arr:np.recarray, contigs_positions_dict:dict):
    arr = arr[np.argsort(arr["position"])]
    float_features = ['norm_mean', 'norm_std', 'dwell_time']
    float_dtypes = [('norm_mean', '<f8'), ('norm_std', '<f8'), ('dwell_time', '<f8')]

    float_arr = arr[float_features].astype(float_dtypes).view('<f8').reshape(-1, 3)
    kmer_arr = arr["reference_kmer"].reshape(-1, 1)
    model_kmer_arr = arr["model_kmer"].reshape(-1, 1)
    tx_pos_arr = arr["position"]
    tx_strand_arr = arr["strand"]
    tx_id_arr = arr["contig"]
    read_indices = arr["read_index"]

    #filter events based on loactions of A base in references
    if tx_strand_arr[0] == "+":
        seg_indices = partions_continuous(tx_pos_arr, contigs_positions_dict["A"])
        ref_11mer_arr = np.array(concensus_kmer(kmer_arr, seg_indices))
        ref_float_arr = np.array([float_arr[i[0]:i[1], :] for i in seg_indices]).reshape(-1, 6 * 3)
    else:
        seg_indices = partions_continuous(tx_pos_arr, contigs_positions_dict["T"])
        ref_11mer_arr = np.array(concensus_model_kmer(model_kmer_arr, seg_indices))
        ref_float_arr = np.array(convert_event(float_arr,seg_indices)).reshape(-1, 6 * 3)

    ref_pos_arr = np.array([tx_pos_arr[(i[-1]-1)] for i in seg_indices])
    ref_strand_arr = np.array([tx_strand_arr[(i[-1]-1)] for i in seg_indices])
    ref_id_arr = np.array([tx_id_arr[(i[-1]-1)] for i in seg_indices])
    ref_read_indices = np.array([read_indices[(i[-1]-1)] for i in seg_indices])

    return tuple((ref_id_arr,ref_read_indices,ref_pos_arr,ref_strand_arr,ref_11mer_arr,ref_float_arr))

