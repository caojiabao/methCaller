import os,re
import multiprocessing
import numpy as np
import pandas as pd
from io import StringIO
from .helper import Consumer,end_queue
from .Utils import *
import time

def WeightedAvg(events_str: str) -> np.recarray:
    r"""
    Function to aggregate features from eventalign.tsv on a single contig id and extract mean current level, dwelling time,
    and current standard deviation from each position in each read
    :param events_str: String corresponding to portion within eventalign.tsv to be processed that can be read as pd.DataFrame object
    :return: np_events (np.recarray): A NumPy record array object that contains the extracted features
    """

    f_string = StringIO(events_str)
    eventalign_result = pd.read_csv(f_string,delimiter='\t',
                                    names=['contig','position','reference_kmer','read_index','strand','event_index',
                                           'event_level_mean','event_stdv','event_length','model_kmer','model_mean',
                                           'model_stdv','standardized_level','start_idx','end_idx'])
    f_string.close()

    cond_successfully_eventaligned = eventalign_result['model_kmer'] != "NNNNNN"

    if cond_successfully_eventaligned.sum() != 0:

        eventalign_result = eventalign_result[cond_successfully_eventaligned]

        l = eventalign_result['event_index'].tolist()
        length = len(l)
        flag = []
        for i in range(0, length - 1):
            if l[i] < l[i + 1]:
                flag.append(1)
            else:
                if l[i] > l[i + 1]:
                    flag.append(0)

        flag = list(set(flag))
        if len(flag) != 1:
            print('WRONG event index')
        else:
            if flag[0] == 1:
                eventalign_result.loc[:, 'strand'] = "+"
            else:
                eventalign_result.loc[:, 'strand'] = "-"
        
        eventalign_result.loc[:, 'length'] = pd.to_numeric(eventalign_result['end_idx']) - \
                pd.to_numeric(eventalign_result['start_idx'])
        eventalign_result.loc[:, 'sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) \
                * eventalign_result['length']
        eventalign_result.loc[:, 'sum_norm_std'] = pd.to_numeric(eventalign_result['event_stdv']) \
                * eventalign_result['length']
        eventalign_result.loc[:, 'sum_dwell_time'] = pd.to_numeric(eventalign_result['event_length']) \
                * eventalign_result['length']

        keys = ['read_index', 'contig', 'position', 'strand', 'reference_kmer', 'model_kmer']  # for groupby
        eventalign_result = eventalign_result.groupby(keys)
        sum_norm_mean = eventalign_result['sum_norm_mean'].sum()
        sum_norm_std = eventalign_result["sum_norm_std"].sum()
        sum_dwell_time = eventalign_result["sum_dwell_time"].sum()

        start_idx = eventalign_result['start_idx'].min()
        end_idx = eventalign_result['end_idx'].max()
        total_length = eventalign_result['length'].sum()

        eventalign_result = pd.concat([start_idx,end_idx],axis=1)
        eventalign_result['norm_mean'] = (sum_norm_mean/total_length).round(1)
        eventalign_result["norm_std"] = sum_norm_std / total_length
        eventalign_result["dwell_time"] = sum_dwell_time / total_length
        eventalign_result = eventalign_result.reset_index()

        features = ['contig', 'read_index',
                    'position', 'strand', 'reference_kmer', 'model_kmer',
                    'norm_mean', 'norm_std', 'dwell_time']
        df_events = eventalign_result[features]
        np_events = np.rec.fromrecords(df_events, names=[*df_events])
        return np_events

    return np.array([])

def parallel_preprocess_tx(eventalign_filepath: str, fasta_filepath:str, out_dir:str, n_processes:int, readcount_min: int):
    r"""
    Function to aggregate segments from the same position within individual read and extract mean current level, dwelling time,
    and current standard deviation from each segment and its neighboring 5 segments from each in eventalign.tsv that
    passes certain count requirements as specified by the user for the purpose of training
    :param eventalign_filepath: String filepath to the eventalign.tsv file
    :param fasta_filepath: String filepath to extract the location of A base in genome
    :param out_dir:  String filepath to the output directory of the indexing function
    :param n_processes: Number of processes used for indexing
    :param readcount_min: Minimum required number of reads expressed by each contig to be considered for feature extraction
    :param readcount_max: Maximum reads expressed by a contig that will be processed for feature extraction
    :return: none
    """

    # Create output paths and locks.
    out_paths, locks = dict(),dict()
    for out_filetype in ['info','log']:
        out_paths[out_filetype] = os.path.join(out_dir, 'data.%s' % out_filetype)
        locks[out_filetype] = multiprocessing.Lock()

    # Writing the starting of the files.
    with open(out_paths['info'], 'w', encoding='utf-8') as f:
        f.write('contig,read_index,contig_position,strand,ref_11mer_sequence,features\n') # header

    open(out_paths['log'], 'w', encoding='utf-8').close()

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.
    consumers = [Consumer(task_queue=task_queue, task_function=preprocess_tx, locks=locks) for i in range(n_processes)]

    for process in consumers:
        process.start()

    df_eventalign_index = pd.read_csv(os.path.join(out_dir,'eventalign.index'))
    df_eventalign_index['contig'] = df_eventalign_index['contig']
    tx_ids = df_eventalign_index['contig'].values.tolist()
    tx_ids = list(dict.fromkeys(tx_ids))
    df_eventalign_index = df_eventalign_index.set_index('contig')
    contig_positions = fasta_to_position(fasta_filepath)
    with open(eventalign_filepath, 'r', encoding='utf-8') as eventalign_result:
        for tx_id in tx_ids:
            contig_position = contig_positions[tx_id]
            data_dict = dict()
            readcount = 0
            for _,row in df_eventalign_index.loc[[tx_id]].iterrows():
                read_index,pos_start,pos_end = row['read_index'],row['pos_start'],row['pos_end']
                eventalign_result.seek(pos_start,0)
                events_str = eventalign_result.read(pos_end-pos_start)
                data = WeightedAvg(events_str)
                if data.size > 1:
                    data_dict[read_index] = data
                readcount += 1
                if readcount % 1000 == 0:
                    task_queue.put((tx_id, contig_position, data_dict, out_paths))
                    readcount = 0
                    data_dict = dict()
            if readcount >= readcount_min:
                task_queue.put((tx_id,contig_position,data_dict, out_paths))

    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue, n_processes)

    # Wait for all the tasks to finish.
    task_queue.join()

def preprocess_tx(tx_id: str, contig_position: dict, data_dict: dict, out_paths: dict, locks: dict):
    r"""

    :param tx_id:  contig id of the portion of eventalign.tsv file to be processed
    :param contig_position: Dictionary containing positions of A base in reference
    :param data_dict: Dictionary containing events for each read
    :param out_paths: A dictionary containing filepath for all the output files produced by the index function
    :param locks: A lock object from multiprocessing library that ensures only one process write to the output file at any given time
    :return: none
    """
    s_time = time.time()
    if len(data_dict) == 0:
        return

    features_arrays = []
    reference_kmer_arrays = []
    read_ids = []
    reference_pos_arrays = []
    reference_strand_arrays = []
    contigs_positions_dict = contig_position #dict:{A:[1,3],T:[1,5]}

    for _,events_per_read in data_dict.items():
        events_per_read = filter_events(events_per_read, contigs_positions_dict)
        features_arrays.append(events_per_read[5])
        reference_kmer_arrays.append(events_per_read[4])
        reference_pos_arrays.append(events_per_read[2])
        reference_strand_arrays.append(events_per_read[3])
        read_ids.append(events_per_read[1])

    if len(features_arrays) == 0:
        return

    features_arrays = np.concatenate(features_arrays)
    reference_kmer_arrays = np.concatenate(reference_kmer_arrays)
    reference_pos_arrays = np.concatenate(reference_pos_arrays)
    reference_strand_arrays = np.concatenate(reference_strand_arrays)
    read_ids = np.concatenate(read_ids)

    assert(len(features_arrays) == len(reference_kmer_arrays) == \
            len(reference_pos_arrays) == len(read_ids) == len(reference_strand_arrays))

    end_time = time.time()
    # write to file.
    log_str = '%s: Data preparation ... Done.' %(tx_id)
    process_reads = '%s reads have been processed. Using time: %s secs' %(len(list(set(read_ids.tolist()))), end_time-s_time)
    with locks['info'], open(out_paths['info'],'a', encoding='utf-8') as f:
        for read_id, pos, strand, kmer, feature in zip(read_ids, reference_pos_arrays, reference_strand_arrays,reference_kmer_arrays,features_arrays):
            features = " ".join(list(map(str,feature.tolist())))
            f.write('%s,%d,%d,%s,%s,%s\n' %(tx_id,read_id,pos,strand,kmer,features))

    with locks['log'], open(out_paths['log'],'a', encoding='utf-8') as f:
        f.write(log_str + '\n'+ process_reads + '\n')