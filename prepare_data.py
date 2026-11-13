# prepare_data.py
# Cuts two FASTA files into 31-mers, one per line: the text input that the
# clustering screen (module 2) of the interface reads.
def extract_kmers(fasta_path, output_path, k=31):
    with open(fasta_path, 'r') as f:
        # Skip the first line (the FASTA header) and join the DNA
        dna = "".join([line.strip() for line in f.readlines()[1:]])

    with open(output_path, 'w') as out:
        # Sliding window: cut the DNA into pieces of 31 letters
        for i in range(len(dna) - k + 1):
            kmer = dna[i:i+k]
            # The Rust side filters out 'N', so they can be skipped right here
            if 'N' not in kmer:
                out.write(kmer + '\n')

print("Slicing real genomes...")
# Assuming the FASTA files have been downloaded
extract_kmers("ecoli.fasta", "sample_A.txt")
extract_kmers("lambda.fasta", "sample_B.txt")
print("Files ready for the interface!")
