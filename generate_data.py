import random

# 1. Generate a "healthy genome" of 10,000 letters
base_dna = "".join(random.choices(['A', 'C', 'G', 'T'], k=10000))

with open("healthy_reference.fasta", "w") as f:
    f.write(">Perfect_Reference_Genome\n")
    f.write(base_dna + "\n")

# 2. Copy the genome, but plant a MUTATION in the middle
mutated_dna = base_dna[:5000] + "AAAAA" + base_dna[5005:]

# Save it as a simulated FASTQ (with the top quality, 'I')
with open("patient_sample.fastq", "w") as f:
    f.write("@Patient_Sample_001\n")
    f.write(mutated_dna + "\n")
    f.write("+\n")
    f.write("I" * len(mutated_dna) + "\n")

print("Test files generated successfully!")
