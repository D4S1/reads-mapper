import mmh3
import numpy as np
from typing import List, Tuple, Dict, Set
from datetime import date
import utils
import math

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
def preprocess_w_set(seq: str, wind_size: int, k: int, hash: int, genome: bool = False) -> Set:
    """
    Generate a set of minimizers from the first `wind_size` characters of a sequence.
    """
    window_kmers = string_to_hashed_kmers(seq[:wind_size], k, hash, genome)
    minimizers = set()
    start_pointer = 0

    for wind_start in range(1, len(seq) - wind_size + 1):
        new_kmer = seq[wind_start + wind_size - k: wind_start + wind_size]
        # Update window k-mers
        if genome:
            window_kmers[start_pointer] = (mmh3.hash(new_kmer, hash), wind_start + wind_size - k + 1)
            minimizers.add(min(window_kmers, key=lambda x: x[0]))
        else:
            window_kmers[start_pointer] = mmh3.hash(new_kmer, hash)
            minimizers.add(min(window_kmers))
        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

    return minimizers

def H_map(w_genom: Set[Tuple[int, int]]) -> Dict[int, List[int]]:
    """
    Convert a set of (hash, position) tuples into a dictionary where each hash maps to a list of positions.
    """
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

def mapping_stage_1(read: str, wind_size: int, k: int, hash: int, H: Dict[int, List[int]], m: int) -> List[Tuple[int, int]]:
    """
    Identify genome position ranges for which the first filtering condition is satisfied.
    """
    T = []
    L = []
    for minimizer in preprocess_w_set(read, wind_size, k, hash):
        L.extend(H.get(minimizer, []))
    L.sort()
    print(L)
    last_added = None
    for i in range(len(L) - m + 1):
        j = i + m - 1
        if L[j] - L[i] < len(read):
            if last_added is None or L[i] - last_added != 1:
                start_pos = max(0, L[j] - len(read) + 1)
                T.append((start_pos, L[i]))
            else:
                T[-1] = (T[-1][0], L[i])
            last_added = L[i]
    return T

def mapping_stage_2(read: str, wind_size: int, k: int, hash: int, T: List[Tuple[int, int]], tau: float) -> List[Tuple[int, float]]:
    """
    Refine genome position ranges by applying the second filtering condition.
    """
    P = []
    L_read = preprocess_w_set(read, wind_size, k, hash)
    for x, y in T:
        L_range = preprocess_w_set(genome[x:x+len(read)], wind_size, k, hash)
        JI = len(L_read.intersection(L_range)) / len(L_read.union(L_range))
        if JI >= tau:
            P.append((x, JI))
    return P

if __name__ == "__main__":
    with open('app.log', 'a') as log_f:
        log_f.write(f"{'='*5} Run date: {date.today()} {'='*5}\n")

    genome = next(iter(utils.read_fasta('data/reference20M.fasta').values()))
    read = 'AGTAACTCTGATAGAAATATGCTAGAGAATATAGTGGGAAAATAAACAGTACTGGTGTTT'

    wind_size = 40
    k = 9
    hash = 1234567
    err_max = 0.1
    delta = 0.1
    s = 4
    tau = get_tau(err_max, k, delta)
    m = math.ceil(s * tau)

    M = preprocess_w_set(genome[:6000], wind_size, k, hash, genome=True)
    H = H_map(M)

    T = mapping_stage_1(read, wind_size, k, hash, H, m)
    print(f"Stage 1 Mapping Results: {T}")

    # P = mapping_stage_2(read, wind_size, k, hash, T, tau)
    # print(f"Stage 2 Mapping Results: {P}")
