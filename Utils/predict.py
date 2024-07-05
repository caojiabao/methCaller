import os
import multiprocessing
import numpy as np
import pandas as pd
from .helper import Consumer,end_queue
from argparse import ArgumentParser
#from sklearn.preprocessing import StandardScaler, MinMaxScaler
from Utils.model import device
import torch
import torch.nn.functional as F

def generate_vector(base):
    if base == 'A':
        return [1, 0, 0, 0]
    elif base == 'C':
        return [0, 1, 0, 0]
    elif base == 'T':
        return [0, 0, 1, 0]
    elif base == 'G':
        return [0, 0, 0, 1]
    else:
        return [0, 0, 0, 0]  # Placeholder for other cases

def transform_sequence(sequence):
    kmer_length = 6
    vectors = []

    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i + kmer_length]
        vector = []
        for base in kmer:
            vector.extend(generate_vector(base))
        vectors.append(vector)

    return np.array(vectors).flatten()


def process_predict(batch_data:pd.DataFrame, model, out_paths:dict, args: ArgumentParser, locks: dict):

    ref_11mer = batch_data[["read_index","ref_11mer_sequence"]]
    vectors = ref_11mer['ref_11mer_sequence'].apply(transform_sequence)
    vectors = np.vstack(vectors)
    seq_features = torch.from_numpy(vectors).type(torch.float)
    seq_features = torch.reshape(seq_features, (len(seq_features), 6, 24))

    features_arr = batch_data.iloc[:, 5:].astype(float).values
    signal_features = torch.from_numpy(features_arr).type(torch.float32)
    signal_features = torch.reshape(signal_features, (-1, 6, 3))
    inputs = torch.cat((signal_features,seq_features),dim=2).to(device)

    model.eval()
    with locks['per_read_methylation'], open(out_paths['per_read_methylation'], 'a', encoding='utf-8') as f:
        with torch.no_grad():
            outputs = model(inputs)

            probabilities = F.softmax(outputs.data, dim=1).cpu().numpy()
            df = batch_data.iloc[:, :5]
            df["score"] = probabilities[:, 1]

            _, predicted_labels = torch.max(outputs.data, dim=1)
            predicted_labels = predicted_labels.cpu().numpy()
            # temp = outputs.numpy() > args.read_proba_threshold
            # df["prediction"] = np.where(temp, 1, 0)
            df['prediction'] = predicted_labels

            for line in df.values.tolist():
                contig, read_index, contig_position, strand, ref_11mer_sequence, score, prediction = line
                f.write("%s,%d,%s,%s,%s,%s,%d\n" % (contig, read_index, contig_position, strand, ref_11mer_sequence, score, prediction))

def process_row(row):
    fifth_column = row[5]
    split_values = fifth_column.split()
    return pd.Series(split_values)

def calculate_function(col):
    if col.name in range(0, 18, 3):
        return (col - 45) / 100
    elif col.name in range(1, 18, 3):
        return (col - 0) / 50
    elif col.name in range(2, 18, 3):
        return col / 0.0075
    else:
        return col

def run_predict(data:str, out_dir:str, args:ArgumentParser, model):

    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['per_read_methylation']:
        out_paths[out_filetype] = os.path.join(out_dir, 'results.%s' % out_filetype)
        locks[out_filetype] = multiprocessing.Lock()

    with open(out_paths['per_read_methylation'], 'w', encoding='utf-8') as f:
        f.write('contig,read_index,contig_position,strand,ref_11mer_sequence,methylation_score,methylation_prediction\n')  # header

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=args.n_processes * 2)

    # Create and start consumers.p
    consumers = [Consumer(task_queue=task_queue, task_function=process_predict, locks=locks) for i in
                 range(args.n_processes)]
    for process in consumers:
        process.start()

    ## Load tasks into task_queue. A task is feature information of one read.
    for chunk in pd.read_csv(data, chunksize=args.chunk_size):
        if len(chunk) == 0:
            return

        df = chunk.apply(process_row, axis = 1)
        chunk = pd.concat([chunk.iloc[:, :5], df], axis=1)
        chunk.iloc[:, 5:] = chunk.iloc[:, 5:].astype(float)
        chunk = chunk.apply(calculate_function)
        task_queue.put((chunk, model, out_paths, args))


    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue, args.n_processes)

    # Wait for all the tasks to finish.
    task_queue.join()
