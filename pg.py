from typing import List, Tuple
import utils

def new_mapping_stage_2(w_read: set, read_length: int, T: List[tuple], M: List[Tuple[int, int]], s: int, tau: float) -> List[Tuple[int, float]]:
    """
    Return a list of tuples of genome position ranges for which the second filtering condition is satisfied
    """
    P = []
    for x, y in T:
        i = x
        j = x + read_length
        w_bi = w_genome_i(i, j, M)

        JI = solve_jackard(w_read, w_bi, s)

        if JI >= tau:
            P.append((i, JI))

        # Slide the window across the range and update L incrementally
        while i < y:
            # Remove minimizers from the left of the window
            w_bi.discard(*w_genome_i(i, i + 1, M))
            # Add minimizers from the right of the window
            w_bi.add(*w_genome_i(j, j + 1, M))

            # Calculate Jaccard index for the updated window
            JI = solve_jackard(w_read, w_bi, s)
            if JI >= tau:
                P.append((i + 1, JI))
            
            # Slide the window
            i += 1
            j += 1

    return utils.merge_ranges(P, read_length)

def w_genome_i(p: int, q: int, M: List[Tuple[int, int]]) -> set[int]:
    return set([M[i][0] for i in range(p, q)])

def solve_jackard(w_read: set, w_genome_i: set, s: int) -> float:
    return  len(sketch(w_read.union(w_genome_i), s).intersection(sketch(w_read, s)).intersection(sketch(w_genome_i, s))) / s

def sketch(w_set: set, s: int) -> set:
    """
    Return the set of the s smallest minimizers.
    """
    return set(sorted(w_set)[:s])