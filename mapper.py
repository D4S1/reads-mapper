from fm_index import FmIndex
import time

start_time = time.time()

genome_path = 'data/reference20M.fasta'

with open(genome_path, 'r') as file:
    lines = file.readlines()
genome = "".join(line.strip() for line in lines[1:])

fm_index = FmIndex(genome)
print(fm_index)
print(f'Time: {(time.time() - start_time)/60} min')