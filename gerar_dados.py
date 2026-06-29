import random

# 1. Vamos gerar um "Genoma Saudável" com 10.000 letras
dna_base = "".join(random.choices(['A', 'C', 'G', 'T'], k=10000))

with open("referencia_saudavel.fasta", "w") as f:
    f.write(">Genoma_Referencia_Perfeito\n")
    f.write(dna_base + "\n")

# 2. Vamos copiar o genoma, mas introduzir uma MUTAÇÃO (SNP) no meio
dna_mutado = dna_base[:5000] + "AAAAA" + dna_base[5005:] 

# Salvar como um FASTQ simulado (com qualidade máxima 'I')
with open("amostra_paciente.fastq", "w") as f:
    f.write("@Amostra_Paciente_001\n")
    f.write(dna_mutado + "\n")
    f.write("+\n")
    f.write("I" * len(dna_mutado) + "\n")

print("Ficheiros de teste gerados com sucesso!")