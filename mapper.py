import mmh3
import numpy as np
from typing import List, Tuple, Dict, Set
from datetime import date
import utils
import math
from math import ceil

# Parameters:
# wind_size: window size (~800)
# k: k-mer length (9)
# err_max: maximum error rate
# n_hash: number of hash functions
# delta: delta for Jaccard similarity estimation

# Create hashed k-mers
def string_to_hashed_kmers(seq: str, k: int, hash: int, genome: bool = False) -> List:
    """
    Generate hashed k-mers from a sequence.
    If genome=True, returns a list of (hash, position) tuples; otherwise, just hashes.
    """
    if genome:
        return [(mmh3.hash(seq[i:i+k], hash), i) for i in range(len(seq) - k + 1)]
    return [mmh3.hash(seq[i:i+k], hash) for i in range(len(seq) - k + 1)]

@utils.timed(1)
def w_set_read(seq: str, wind_size: int, k: int, hash: int) -> Set[int]:
    """
    Generate a set of minimizers, corresponding to each `wind_size` characters of a read.
    """
    window_kmers = string_to_hashed_kmers(seq[:wind_size], k, hash)
    minimizers = set()
    start_pointer = 0

    for wind_start in range(1, len(seq) - wind_size + 1):
        new_kmer = seq[wind_start + wind_size - k: wind_start + wind_size]
        # Update window k-mers
        window_kmers[start_pointer] = mmh3.hash(new_kmer, hash)
        minimizers.add(min(window_kmers))
        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

    return minimizers

@utils.timed(1)
def w_set_genome(seq: str, wind_size: int, k: int, hash: int) -> List[Tuple[int, int]]:
    """
    Generate a list of minimizers, along with their positions, corresponding to each `wind_size` characters of the genome.
    """
    window_kmers = string_to_hashed_kmers(seq[:wind_size], k, hash, genome=True)
    minimizers = []
    start_pointer = 0

    for wind_start in range(1, len(seq) - wind_size + 1):
        new_kmer = seq[wind_start + wind_size - k: wind_start + wind_size]
        # Update window k-mers
        window_kmers[start_pointer] = (mmh3.hash(new_kmer, hash), wind_start + wind_size - k + 1)
        minimizers.append(min(window_kmers, key=lambda x: x[0]))
        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

    return minimizers

def H_map(w_genom: List[Tuple[int, int]]) -> Dict[int, List[int]]:
    """
    Convert a set of (hash, position) tuples into a dictionary where each hash maps to a list of positions.
    """
    w_genom = set(w_genom)
    result_dict = {}
    for h, pos in w_genom:
        result_dict.setdefault(h, []).append(pos)
    return result_dict

def sketch(w_set: Set, s: int) -> Set:
    """
    Return the set of the s smallest minimizers.
    """
    return set(sorted(w_set)[:s])

def winnowed_minhash_estimate(w_read: Set, w_genome_i: Set, s: int) -> float:
    """
    Estimate Jaccard similarity using winnowed minhash.
    """
    s_w_genome_i = sketch(w_genome_i, s)
    s_w_read = sketch(w_read, s)
    s_union = sketch(w_read.union(w_genome_i), s)
    return len(s_union.intersection(s_w_read).intersection(s_w_genome_i)) / len(s_union)

def get_tau(err_max: float, k: int, delta: float) -> float:
    """
    Calculate Jaccard similarity threshold given error rate and other parameters.
    """
    return 1 / (2 * math.exp(err_max * k) - 1) - delta

def mapping_stage_1(w_read: set, read_length: int, H: Dict[int, List[int]], m: int) -> List[Tuple[int, int]]:
    """
    Identify genome position ranges for which the first filtering condition is satisfied.
    """
    T = []
    L = []
    for minimizer in w_read:
        L.extend(H.get(minimizer, []))
    L.sort()
    last_added = None
    for i in range(len(L) - m + 1):
        j = i + m - 1
        if L[j] - L[i] < read_length:
            if last_added is None or L[i] - last_added != 1:
                start_pos = max(0, L[j] - read_length + 1)
                T.append((start_pos, L[i]))
            else:
                T[-1] = (T[-1][0], L[i])
            last_added = L[i]
    return T

def mapping_stage_2(w_read: set, read_length: int, T: List[tuple], M: List[Tuple[int, int]], s: int, tau: float) -> List[Tuple[int, float]]:
    """
    Returns a list of tuples of genome position ranges
    for which the second filtering condition is satisfied
    """
    P = []
    L = {h: 1 for h in w_read}  # initialize the map L with read minimizers
    for x, y in T:
        i = x
        j = x + read_length
        w_bi = get_minimizers(i, j, M)
        for h in w_bi:
            if h in L:
                L[h] = 1  # shared minimizer between read and Bi
            else:
                L[h] = 0  # unique to Bi
        JI = solve_jackard(L, s)
        if JI >= tau:
            P.append((i, JI))

        # Slide the window across the range and update L incrementally
        while i < y:
            # Remove minimizers from the left of the window
            left_minimizers = get_minimizers(i, i + 1, M)
            for h in left_minimizers:
                if h in L:
                    del L[h]

            # Add minimizers from the right of the window
            right_minimizers = get_minimizers(j, j + 1, M)
            for h in right_minimizers:
                if h in w_read:
                    L[h] = 1  # Shared minimizer
                else:
                    L[h] = 0  # Unique minimizer
            
            # Calculate Jaccard index for the updated window
            JI = solve_jackard(L, s)
            if JI >= tau:
                P.append((i + 1, JI))
            
            # Slide the window
            i += 1
            j += 1

    return P

def get_minimizers(p: int, q: int, M: List[Tuple[int, int]]) -> Set[int]:
    return set([M[i][0] for i in range(p, q)])

def solve_jackard(L: dict, s: int) -> float:
    shared = sum(v for v in L.values())
    return min(shared, s) / s

def mapper(read: str, M: list, H: set, wind_size: int, k: int, hash: int, err_max: float, delta:float) -> List:

    w_read = w_set_read(read, wind_size, k, hash)
    s = len(w_read)
    tau = get_tau(err_max, k, delta)
    m = ceil(s * tau)
    read_length = len(read)

    T = mapping_stage_1(w_read, read_length, H, m)
    P = mapping_stage_2(w_read, read_length, T, M, s, tau)

    return P


def main(reads_filename, genome_filename, wind_size, k, hash, err_max, delta):

    reads = utils.read_fasta(reads_filename)
    genome = next(iter(utils.read_fasta(genome_filename).values()))

    M = w_set_genome(genome[:6000], wind_size, k, hash)
    H = H_map(M)

    for id, read in reads.items():
        P = mapper(read, M, H, wind_size, k, hash, err_max, delta)
        print(f'{id=}\n{P}')

if __name__ == "__main__":
    with open('app.log', 'a') as log_f:
        log_f.write(f"{'='*5} Run date: {date.today()} {'='*5}\n")

    wind_size = 40
    k = 9
    hash = 1234567
    err_max = 0.1
    delta = 0.1

    main('data/reads_test.fasta', 'data/reference20M.fasta', wind_size, k, hash, err_max, delta)
