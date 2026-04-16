import concurrent.futures
import os

fasta_file = "dereplicated_contigs.fasta" 
id_file = "viral_list.txt"  
output_file = "vOTUs.fasta" 
num_threads = 30


with open(id_file, "r") as f:
    id_list = {line.strip() for line in f}


def process_chunk(lines, id_list):
    sequences = {}
    current_id = ""
    for line in lines:
        if line.startswith(">"):
            current_id = line[1:].strip()
            if current_id in id_list:
                sequences[current_id] = ""
        else:
            if current_id in id_list:
                sequences[current_id] += line.strip()
    return sequences


chunk_size = 1000 
with open(fasta_file, "r") as f:
    lines = f.readlines()

chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]


all_sequences = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = [executor.submit(process_chunk, chunk, id_list) for chunk in chunks]
    for future in concurrent.futures.as_completed(futures):
        sequences = future.result()
        all_sequences.update(sequences)


with open(output_file, "w") as f:
    for seq_id, sequence in all_sequences.items():
        f.write(">" + seq_id + "\n" + sequence + "\n")
