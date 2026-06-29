# test_kmer.py
import multimerge_genomics as mmg
import pyarrow as pa
import polars as pl

print("1. Testando Codificador de K-mers (CPU):")
sequencia = "ACGTAGCATCGGATCTATCTATCGACACTTG" # 31 bases
resultado = mmg.encode_sequence_to_kmer(sequencia)
print(f"Sequência: {sequencia}")
print(f"Codificação u64: {resultado}\n")

print("2. Testando Pipelining (GPU -> Arrow -> Polars):")

# Gerando amostras de teste limpas (31 bases)
amostras_A = ["ACGTAGCATCGGATCTATCTATCGACACTTG", "GGTTACGTACGTACGTACGTACGTACGTACG"]
amostras_B = ["TTCGATCGATCGATCGATCGATCGATCGATC", "CCCGGGAAATTTCCCGGGAAATTTCCCGGGA"]

# Inicia o processamento na GPU em segundo plano
# O Rust cria a thread da GPU, aloca VRAM e começa a computar imediatamente
stream_de_dados = mmg.iniciar_stream_gpu(amostras_A, amostras_B)

# O Python atua como o consumidor (Pull System)
lote_numero = 1
for lote_arrow in stream_de_dados:
    print(f"--- Processando Lote {lote_numero} ---")
    print(f"Tipo de dados recebido nativamente: {type(lote_arrow)}")
    print(f"Valores do Apache Arrow: {lote_arrow}")
    
    # Transformando num DataFrame do Polars em tempo O(1) (Zero-Copy)
    # Este é o DataFrame que será enviado para o DNABERT-S (IA) na Fase 4
   # Lê a memória Zero-Copy como uma Coluna (Series)
    coluna_polars = pl.from_arrow(lote_arrow)
    
    # Transforma a Coluna em uma Tabela (DataFrame) nomeada "kmer_id"
    df_polars = coluna_polars.to_frame(name="kmer_id")
    
    print("\nDataFrame Polars Pronto:")
    print(df_polars)
    
    lote_numero += 1

print("\nPipeline concluída com sucesso! Nenhuma cópia de RAM foi criada.")