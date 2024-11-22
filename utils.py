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
