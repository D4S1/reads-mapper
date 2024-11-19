
import mmh3
import numpy as np
from typing import List

# parametry:

# wind_size: windows size (~800)
# k: kmer length (9)
# err: error rate 
# n_hash: number of has functions
# dt_est: delta for JC estimation


# create hashed kmers
def string_to_hashed_kmers(seq: str, k: int, n_hash: int) -> np.array:
    # matrix shape (wind_size - k +1, n_hash)
    return np.array([hash_kmer(seq[i:i+k], n_hash) for i in range(len(seq)-k+1)])

def hash_kmer(kmer: str, n_hash: int) -> List:
    return [mmh3.hash(kmer, seed) for seed in range(n_hash)]
    

def preprocess_genome(genome, wind_size, k, n_hash):

    window_kmers = string_to_hashed_kmers(genome[:wind_size], k, n_hash)
    print(f'{window_kmers.shape=}')
    minimazers = [set(window_kmers.min(axis=0))]
    start_pointer = 0
    for wind_start in range(1, len(genome) - wind_size):
        
        # extract kmer that appear by 
        new_kmer = genome[wind_start + wind_size - k: wind_start + wind_size]

        # update a window kmers
        window_kmers[start_pointer] = hash_kmer(new_kmer, n_hash)
        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

        minimazers.append(set(window_kmers.min(axis=0)))

    return minimazers



seq = "AAAAACCCCCCGGGGTTTTAAA"
wind_size = 10
k = 5
n_hash = 10

# mat = string_to_hashed_kmers(seq=seq[:10], k=k, n_hash=n_hash)

print(preprocess_genome(seq, 10, 5, 3))
