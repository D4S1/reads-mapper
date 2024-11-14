from fm_index import FmIndex
from alignment import align_banded_global
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
    mapped_kmers = []

    for ith_kmer in range(101):
        kmer = read[ith_kmer*size: size*(ith_kmer + 1)]

        first, last = ref_genome.range(kmer)

        for o in range(first, last+1):
            real_idx = ref_genome.resolve(o)
            mapped_kmers.append({ith_kmer*size: [kmer, real_idx]})

    return mapped_kmers


def main(reads, ref_genome):
    for id, read in reads.items():
        seeds = map_seeds(read, ref_genome, size=19)
    
    return seeds

dummy_data = {'read_0':'ATCCTATGAAAAATGCAAATGTGTCCATAAGCACAACTCTGTAATTATAGCTCTGCAGGTTAACTGATCTGAGTCTGAGTTTCCTCATTTGAATAACTGGGGTACTAATAACACCTATTTCACATAGATGATGTGAATATTCTATAAGCAATAAGTAGAAAATTAGCTGTTATACATAGGAAGATGAACAAATGATGGCTGTTGTTAAAAAAATAAAAAGCAAAGAGAGTAATTCTTTTTTTAAAAAGCTTTAATGAGGATTAAGTGACCTACTATCTGCAAGAGGCTTTGCACCGAGTAGGAACTCAGCAACATTTATTTCCTGCACTCCCTTTCTTCCTTCTTTGCTTTGTGAAAAAAAGACAAGCTGATGGACACACTCATTAGGGATAATTTATACTATTTTCTGAAGCAGAGTAGACTCAGTCATTTCTTTTATCCAGTTGCTAGTAATAGATGTTTGGCAGAAAAAGTCTCACTCAATGATGCAAGAGGTAACTGACTCACTAAAGACGTGGCATGCAAAATCCCGTTGGTGTCCATCCAGTTACATTTCTCCCCCAATGGTGGATGCCTAGAATTTCACAGTCTCTTCTCTACTTTTGATTTTGGTACCTCTTAAACCATTACTGCTCCTACAAATTCGTCTGAGTAATTATCCATACCTTCCCTGGAGTATGAGAGTGAAATAGATCTATAAGTAGATTTTTAGAATTAATGCTTAGGATTTGCATGAGGGAAGGAAGAGAGGAGGAAAGGAATGCTTGTGGCTGCTGATTTCAGAAGGAAATGGAAGTAGAATCTACAAACAGGGAAGGGAAACATATTGAAGACAAAAAAACAGTTTCCTACTGCATGAATTTTCCTGGATTTAATCCTTAATCCGAATGCATCAGGAACATCTACCATGCATTTGAGATGTCTGAGTGGGTCCAATCCAAGGGAAAGAAGCAAGTACTTTCTTAATGGTTTCTTTGTTGCTACATTAGCATTCCAAACT'}

if __name__ == "__main__":
    # genome_fm_index = preprocess_genome(genome_path)
    with open('preprocessed_genome.pkl', 'rb') as f:
        ref_genome = pickle.load(f)

    res = main(dummy_data, ref_genome)
    print(f'{len(res)=}')
