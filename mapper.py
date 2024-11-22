
import mmh3
import numpy as np
from typing import List
from datetime import date
import utils

# parametry:

# wind_size: windows size (~800)
# k: kmer length (9)
# err: error rate 
# n_hash: number of has functions
# dt_est: delta for JC estimation


# create hashed kmers
def string_to_hashed_kmers(seq: str, k: int, hash: int) -> List:
    # list len: wind_size - k + 1
    return [(mmh3.hash(seq[i:i+k], hash), i) for i in range(len(seq)-k+1)]
    
@utils.timed(1)
def preprocess_genome(genome, wind_size, k, hash):

    window_kmers = string_to_hashed_kmers(genome[:wind_size], k, hash)
    minimazers = set()
    start_pointer = 0

    for wind_start in range(1, len(genome) - wind_size):
        
        # extract kmer that appear by 
        new_kmer = genome[wind_start + wind_size - k: wind_start + wind_size]

        # update a window kmers
        window_kmers[start_pointer] = (mmh3.hash(new_kmer, hash), wind_start + wind_size - k +1)
        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

        minimazers.add(min(window_kmers, key=lambda x: x[0]))

    return minimazers

def sketch(w_set:set) -> set:
    "Returns s smallest minimazers"
    pass

def winnowed_minhash_estimate(w_read: set, w_genome_i:set) -> float:
    "Returns estimation of Jaccard similarity"
    pass

def tau(err_max:float, k:int, delta:float) -> float:
    "Returns Jaccard similarity estimation given error rate"
    pass


if __name__ == "__main__":
    # Run log init
    with open('app.log', 'a') as log_f:
        log_f.write(f"{'='*5}Run date: {date.today()}{'='*5}\n")

    genome = next(iter(utils.read_fasta('data/reference20M.fasta').values()))

    wind_size = 10
    k = 9
    n_hash = 1

    res = preprocess_genome(genome[:20], wind_size, k, n_hash)
    print(res)