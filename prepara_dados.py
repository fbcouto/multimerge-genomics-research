# prepara_dados.py
def extrair_kmers(arquivo_fasta, arquivo_saida, k=31):
    with open(arquivo_fasta, 'r') as f:
        # Ignora a primeira linha (cabeçalho do FASTA) e junta o DNA
        dna = "".join([linha.strip() for linha in f.readlines()[1:]])
    
    with open(arquivo_saida, 'w') as out:
        # Janela deslizante: fatia o DNA em pedaços de 31 letras
        for i in range(len(dna) - k + 1):
            kmer = dna[i:i+k]
            # O seu Rust filtra 'N', então podemos ignorá-los logo aqui
            if 'N' not in kmer:
                out.write(kmer + '\n')

print("Fatiando Genomas reais...")
# Supondo que baixou os FASTAs
extrair_kmers("ecoli.fasta", "amostra_A.txt")
extrair_kmers("lambda.fasta", "amostra_B.txt")
print("Arquivos prontos para o Marimo!")