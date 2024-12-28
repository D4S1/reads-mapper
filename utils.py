import logging
import time
import pickle
import argparse

logging.basicConfig(
    level=logging.DEBUG,
    format=f"%(message)s",
    handlers=[
        logging.FileHandler("app.log"),  # Save logs to a file named 'app.log'
        # logging.StreamHandler()         # Also print logs to the console
    ]
)
logger = logging.getLogger("my-logger")

def timed(num_ops):
    """
    A decorator to log execution time for the decorated function.
    It also calculates average time per operation if `num_operations` > 1.
    """
    def decorator(func):

        def wrapper(*args, **kwargs):
            start = time.time()
            result = func(*args, **kwargs)
            end = time.time()

            total_time = end - start

            # Log total time
            logger.debug(f"{func.__name__} ran in {round(total_time, 2)}s")

            # Log average time per operation if num_ops > 1
            if num_ops > 1:
                avg_time = total_time / num_ops
                logger.debug(
                    f"Time(s) per read: {round(avg_time, 2)}s"
                )
            return result
        return wrapper
    return decorator

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_id = None
        sequence_data = []
        
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = ''.join(sequence_data)
                sequence_id = line[1:]  # Remove the '>' character
                sequence_data = []
            else:
                sequence_data.append(line)
        
        # Add the last sequence to the dictionary
        if sequence_id:
            sequences[sequence_id] = ''.join(sequence_data)
    
    return sequences

def save_to_file(obj, filename):
    with open(filename, 'wb') as file:
        pickle.dump(obj, file)

def load_pickle(filename):
    with open(filename, 'rb') as file:
        return pickle.load(file)
    
def merge_ranges(ranges, read_len):
    window_size = 100
    last_end = ranges[0][0] + read_len + window_size
    reg_n = 1
    
    mereged_ranges = [(ranges[0][0], last_end, ranges[0][1])]

    for current_start, jc in ranges[1:]:
        current_end = current_start + read_len + window_size

        if current_start <= last_end + 1:
            mereged_ranges[-1] = (mereged_ranges[-1][0], current_end, mereged_ranges[-1][2] + jc)
            reg_n += 1

        else:
            mereged_ranges[-1] = (mereged_ranges[-1][0], mereged_ranges[-1][1], mereged_ranges[-1][2] / reg_n)
            reg_n = 1
            mereged_ranges.append((current_start, current_end, jc))

        last_end = current_end
    
    mereged_ranges[-1] = (mereged_ranges[-1][0], mereged_ranges[-1][1], mereged_ranges[-1][2] / reg_n)

    return sorted(mereged_ranges, key=lambda x:-x[2])


def accuracy(out_file, test_file):
    
    with open(out_file, 'r') as file:
        out_locs = [list(map(int, line.split('\t')[1:])) for line in file.read().strip().split('\n')]

    with open(test_file, 'r') as file:
        test_locs = [list(map(int, line.split('\t')[1:])) for line in file.read().strip().split('\n')]

    acc = 0
    mapped = 0
    for id, (pre_loc, test_loc) in enumerate(zip(out_locs, test_locs)):
        if pre_loc == [0, 0]:
            print('skipped')
            continue
        mapped += 1
        if score(pre_loc, test_loc):
            acc += 1
        else:
            print(f'badly mapped read id: {id}')

    return acc / mapped


import pandas as pd

def accuracy(out_file, test_file):
    # Load files into DataFrames
    out_df = pd.read_csv(out_file, sep='\t', header=None, names=['id', 'start', 'end'])
    test_df = pd.read_csv(test_file, sep='\t', header=None, names=['id', 'start', 'end'])

    out_df['id'] = out_df['id'].astype(str).str.strip()
    test_df['id'] = test_df['id'].astype(str).str.strip()

    # Merge DataFrames on the 'id' column
    merged_df = out_df.merge(test_df, on='id', how='left', suffixes=('_pred', '_test'))

    acc = 0
    for _, row in merged_df.iterrows():
        test_loc = [row['start_test'], row['end_test']]
        pre_loc = [row['start_pred'], row['end_pred']]

        if score(pre_loc, test_loc):
            acc += 1
        else:
            print(f'badly mapped read id: {row["id"]}')

    print(f'mapped {len(merged_df)*100/len(test_df)}%')
    return acc / len(merged_df)


def mapped_reads_and_acc(out_file, test_file):

    with open(out_file, 'r') as file:
        out_data = [line.strip().split('\t') for line in file.readlines()]
        out_dict = {line[0]: list(map(int, line[1:])) for line in out_data}  # dict: id -> [start, end]

    with open(test_file, 'r') as file:
        test_data = [line.strip().split('\t') for line in file.readlines()]
        test_dict = {line[0]: list(map(int, line[1:])) for line in test_data}  # dict: id -> [start, end]

    mapped_reads = len(out_dict) / len(test_dict) if test_dict else 0

    acc = 0
    for read_id, pre_loc in out_dict.items():
            test_loc = test_dict[read_id]
            if score(pre_loc, test_loc):
                acc += 1
    acc = acc / len(out_dict) if out_dict else 0.0

    return mapped_reads, acc


def score(pre_loc, test_loc):
    return abs(pre_loc[0] - test_loc[0]) + abs(pre_loc[1] - test_loc[1]) <= 20
