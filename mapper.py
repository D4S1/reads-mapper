from fm_index import FmIndex
from alignment import align_banded
import time
import pickle

start_time = time.time()

genome_path = 'data/reference20M.fasta'
reads_path = 'data/reads20Mb.fasta'

def preprocess_genome(genome_path):
    with open(genome_path, 'r') as file:
        lines = file.readlines()
    genome = "".join(line.strip() for line in lines[1:])

    fm_index = FmIndex(genome)
    print(f'Time: {(time.time() - start_time)/60} min')
    fm_index.save_to_file('preprocessed_genome.pkl')
    return fm_index


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


def map_seeds(read, ref_genome, size=9):
    mapped_kmers = {}

    for ith_kmer in range(101):
        kmer = read[ith_kmer*size: size*(ith_kmer + 1)]
        mapped_kmers[ith_kmer*size] = [kmer, []]

        first, last = ref_genome.range(kmer)
        for o in range(first, last+1):
            real_idx = ref_genome.resolve(o)
            mapped_kmers[ith_kmer*size][1].append(real_idx)

    return mapped_kmers


def main(reads, ref_genome, size=9):
    for id, read in reads.items():
        seeds = map_seeds(read, ref_genome, size=size)

        for read_idx, (seed, ref_indices)  in seeds.items():
            for ref_idx in ref_indices:
                pre_seq, suf_seq = read[:read_idx], read[read_idx+size:]
                pre_ref, suf_ref = ref_genome.seq[ref_idx - 100 - read_idx :ref_idx], ref_genome.seq[ref_idx+size: ref_idx+size + len(suf_seq) + 100]
                res = align_banded(pre_seq,pre_ref, k=100)
                print(res)
    return res

dummy_data = {'read_0':'ATCCTATGAAAAATGCAAATGTGTCCATAAGCACAACTCTGTAATTATAGCTCTGCAGGTTAACTGATCTGAGTCTGAGTTTCCTCATTTGAATAACTGGGGTACTAATAACACCTATTTCACATAGATGATGTGAATATTCTATAAGCAATAAGTAGAAAATTAGCTGTTATACATAGGAAGATGAACAAATGATGGCTGTTGTTAAAAAAATAAAAAGCAAAGAGAGTAATTCTTTTTTTAAAAAGCTTTAATGAGGATTAAGTGACCTACTATCTGCAAGAGGCTTTGCACCGAGTAGGAACTCAGCAACATTTATTTCCTGCACTCCCTTTCTTCCTTCTTTGCTTTGTGAAAAAAAGACAAGCTGATGGACACACTCATTAGGGATAATTTATACTATTTTCTGAAGCAGAGTAGACTCAGTCATTTCTTTTATCCAGTTGCTAGTAATAGATGTTTGGCAGAAAAAGTCTCACTCAATGATGCAAGAGGTAACTGACTCACTAAAGACGTGGCATGCAAAATCCCGTTGGTGTCCATCCAGTTACATTTCTCCCCCAATGGTGGATGCCTAGAATTTCACAGTCTCTTCTCTACTTTTGATTTTGGTACCTCTTAAACCATTACTGCTCCTACAAATTCGTCTGAGTAATTATCCATACCTTCCCTGGAGTATGAGAGTGAAATAGATCTATAAGTAGATTTTTAGAATTAATGCTTAGGATTTGCATGAGGGAAGGAAGAGAGGAGGAAAGGAATGCTTGTGGCTGCTGATTTCAGAAGGAAATGGAAGTAGAATCTACAAACAGGGAAGGGAAACATATTGAAGACAAAAAAACAGTTTCCTACTGCATGAATTTTCCTGGATTTAATCCTTAATCCGAATGCATCAGGAACATCTACCATGCATTTGAGATGTCTGAGTGGGTCCAATCCAAGGGAAAGAAGCAAGTACTTTCTTAATGGTTTCTTTGTTGCTACATTAGCATTCCAAACT'}

if __name__ == "__main__":
    # genome_fm_index = preprocess_genome(genome_path)
    with open('preprocessed_genome.pkl', 'rb') as f:
        ref_genome = pickle.load(f)

    s = time.time()
    res = main(dummy_data, ref_genome)
    print(f'Time: {time.time() - s}')
    print(f'{len(res)=}')
