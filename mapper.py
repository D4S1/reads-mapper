from fm_index import FmIndex
from alignment import align_banded_global
import time

start_time = time.time()

genome_path = 'data/reference20M.fasta'
reads_path = 'data/reads20Mb.fasta'

def preprocess_genome(genome_path):
    with open(genome_path, 'r') as file:
        lines = file.readlines()
    genome = "".join(line.strip() for line in lines[1:])

    fm_index = FmIndex(genome)
    print(f'Time: {(time.time() - start_time)/60} min')
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

genome_fm_index = preprocess_genome(genome_path)
sequences = read_fasta(reads_path)


def map_seed(read, ref_genome, size=9):

    kmers = [{i*size: [read[i*size:size(1+1)], None]} for i in range(101)]

    for kmer in kmers:
        first, _ = ref_genome.range(kmers[kmer][0])
        real_idx = ref_genome.resolve(first)
        kmers[kmer][1] = real_idx

    return kmers


def main(reads_file, ref_file):

    reads = read_fasta(reads_file)
    
    with open(ref_file, 'r') as f:
        ref_genome = f


genome_fm_index = preprocess_genome(genome_path)
print(genome_fm_index)