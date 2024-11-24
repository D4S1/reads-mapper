import logging
import time

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


def merge_ranges(ranges, read_len):

    last_end = ranges[0][0] + read_len
    reg_n = 1
    
    mereged_ranges = [(ranges[0][0], last_end, ranges[0][1])]

    for current_start, jc in ranges[1:]:
        current_end = current_start + read_len

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
        out_locs = [list(map(int, line.split('\t'))) for line in file.read().split('\n')]

    with open(test_file, 'r') as file:
        test_locs = [list(map(int, line.split('\t'))) for line in file.read().split('\n')]

    acc = 0

    for pre_loc, test_loc in zip(out_locs, test_locs):

        if score(pre_loc, test_loc):
            acc += 1
    
    return acc / len(out_locs)



def score(pre_loc, test_loc):
    return abs(pre_loc[0] - test_loc[0]) + abs(pre_loc[1] - test_loc[1]) <= 20
