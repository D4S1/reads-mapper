import mmh3
import numpy as np
from typing import List
from datetime import date
import utils
import math

# parametry:

# wind_size: windows size (~800)
# k: kmer length (9)
# err: error rate 
# n_hash: number of has functions
# dt_est: delta for JC estimation


# create hashed kmers
def string_to_hashed_kmers(seq: str, k: int, hash: int, genome: bool = False) -> List:
    # list len: wind_size - k + 1
    if genome:
        return [(mmh3.hash(seq[i:i+k], hash), i) for i in range(len(seq)-k+1)]
    return[mmh3.hash(seq[i:i+k], hash) for i in range(len(seq)-k+1)]
    
@utils.timed(1)
def preprocess_w_set(seq: str, wind_size: int, k: int, hash: int, genome: bool = False) -> set:

    window_kmers = string_to_hashed_kmers(seq[:wind_size], k, hash, genome)
    minimizers = set()
    start_pointer = 0

    for wind_start in range(1, len(seq) - wind_size + 1):
        
        # extract kmer that appear by 
        new_kmer = seq[wind_start + wind_size - k: wind_start + wind_size]

        # update a window kmers
        if genome: 
            window_kmers[start_pointer] = (mmh3.hash(new_kmer, hash), wind_start + wind_size - k +1)
            minimizers.add(min(window_kmers, key=lambda x: x[0]))
        else: 
            window_kmers[start_pointer] = mmh3.hash(new_kmer, hash)
            minimizers.add(min(window_kmers))

        start_pointer = (start_pointer + 1) % (wind_size - k + 1)

    return minimizers

def H_map(w_genom: set) -> dict:
    """
    Convert a set of (h, pos) tuples into a dictionary where each h key
    maps to a list of its corresponding pos values.
    """
    result_dict = {}
    for h, pos in w_genom:
        # If h is already a key, append pos to its list
        if h in result_dict:
            result_dict[h].append(pos)
        # If h is not a key, add it with pos as the start of its list
        else:
            result_dict[h] = [pos]
    return result_dict

def sketch(w_set: set, s: int) -> set:
    "Returns set of s smallest minimazers"
    return set(sorted(w_set)[:s])

def winnowed_minhash_estimate(w_read: set, w_genome_i: set, s: int) -> float:
    "Returns estimation of Jaccard similarity"
    s_w_genome_i = sketch(w_genome_i, s)
    s_w_read = sketch(w_read, s)
    s_union = sketch(w_read.union(w_genome_i), s)
    return len(s_union.intersection(s_w_read).intersection(s_w_genome_i)) / len(s_union)

def get_tau(err_max: float, k: int, delta: float) -> float:
    "Returns Jaccard similarity estimation given error rate"
    return 1/(2 * math.exp(err_max * k) - 1) - delta

def get_minimizers(i: int, j: int):
    "Returns hash values for positions between i and j in the genome"
    # potrzebujemy tutaj mieć jednak listę minimizerów, zeby łatwo móc się do nich dostać
    pass

def solve_jackard(L: set, s: int) -> float:
    return sum(L)/s  # to jest chyba źle

def mapping_stage_1(read: str, wind_size: int, k: int, hash: int, H: dict, m: int) -> List[tuple]:
    """
    Returns a list of tuples of genome position beginning ranges
    for which the first filtering condition is satisfied
    """

    # to do: mozna zoptymalizować nakładające się range

    T, L = [], []
    for minimizer in preprocess_w_set(read, wind_size, k, hash):
        try:
            L.extend(H[minimizer])
        except KeyError:
            continue
    L.sort()
    print(f'{L=}')
    last_added = None
    for i in range(len(L) - m + 1):
        j = i + m - 1
        if L[j] - L[i] < len(read):
            if last_added is None or L[i] - last_added != 1:
                if L[j] - len(read) + 1 < 0:
                    T.append([0, L[i]])
                else:
                    T.append([L[j] - len(read) + 1, L[i]])
                last_added = L[i]

            else:
                # if the difference between the last added value and current value is one,
                # we can extend the existing range instead of creating a new one
                T[-1][-1] = L[i]
                last_added = L[i]
    return T

def mapping_stage_2(read: str, wind_size: int, k: int, hash: int, T: List[tuple]):
    """
    Returns a list of tuples of genome position ranges
    for which the second filtering condition is satisfied
    """
    L = {}
    L0 = {preprocess_w_set(read, wind_size, k, hash)}
    P = []
    for x, y in T:
        i = x
        j = x + len(read)
        L = L0
        L.insert(get_minimizers(i, j))
        JI = solve_jackard(L)
        if JI >= tau:
            P.append((i, JI))
        while i <= y:
            L.detete(get_minimizers(i, i+1))
            L.insert(get_minimizers(j, j+1))
            JI = solve_jackard(L)
            if JI >= tau:
                P.append((i, JI))
            i += 1
            j += 1
    return P

def filter_step_1(w_read: set, w_genome_i: set, s: int, err_max: float, k: int, delta: float) -> set:
    "Returns the set of potential mapping position after stage 1 filtering"


def filter_step_2(w_read: set, w_genome_i: set, s: int, err_max: float, k: int, delta: float) -> set:
    pass

if __name__ == "__main__":
    # Run log init
    with open('app.log', 'a') as log_f:
        log_f.write(f"{'='*5}Run date: {date.today()}{'='*5}\n")

    genome = next(iter(utils.read_fasta('data/reference20M.fasta').values()))
    # read = next(iter(utils.read_fasta('data/reads20Ma.fasta').values()))[:10]
    read = 'AGTAACTCTGATAGAAATATGCTAGAGAATATAGTGGGAAAATAAACAGTACTGGTGTTT'

    wind_size = 40
    k = 9
    hash = 1234567
    err_max = 0.1
    delta = 0.1
    s = 4
    read_length = len(read)
    tau = 1/(2 * math.exp(err_max * k) - 1) - delta
    m = math.ceil(s * tau)
    # m = 2

    M = preprocess_w_set(genome[:6000], wind_size, k, hash, genome=True)
    H = H_map(M)

    T = mapping_stage_1(read, wind_size, k, hash, H, m)
    print(T)

