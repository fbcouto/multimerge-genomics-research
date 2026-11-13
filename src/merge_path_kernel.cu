#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// =========================================================
// CONTRATO COM O RUST
// =========================================================
// Todo lancador devolve 0 quando deu certo, ou um codigo diferente de 0:
// o codigo de erro do CUDA (positivo) ou MM_ERRO_TAMANHO (negativo) quando a
// entrada nao cabe nos indices do kernel. Quem chama (src/gpu.rs) cai para a
// CPU quando nao e 0.
//
// Antes nada era conferido: sem placa, sem driver, ou com dado maior que a
// VRAM, o cudaMalloc falhava calado e o resultado era lixo ou queda do
// processo.
#define MM_ERRO_TAMANHO (-1)
#define MM_TENTA(chamada) do { cudaError_t mm_e = (chamada); if (mm_e != cudaSuccess) { erro = (int)mm_e; goto fim; } } while (0)

extern "C" {
    // Quantas placas CUDA este processo enxerga (0 sem placa ou sem driver).
    int mm_cuda_dispositivos() {
        int n = 0;
        if (cudaGetDeviceCount(&n) != cudaSuccess) {
            cudaGetLastError(); // limpa o erro pendente
            return 0;
        }
        return n;
    }

    // Memoria livre e total da placa em uso, em bytes. Devolve 0 se deu certo.
    int mm_cuda_memoria(size_t* livre, size_t* total) {
        cudaError_t e = cudaMemGetInfo(livre, total);
        if (e != cudaSuccess) {
            cudaGetLastError();
            return (int)e;
        }
        return 0;
    }
}

// 1. Busca Binária do Caminho de Fusão
__device__ int merge_path_search(const uint64_t* A, int lenA, const uint64_t* B, int lenB, int diag) {
    int begin = max(0, diag - lenB);
    int end = min(diag, lenA);
    while (begin < end) {
        int mid = (begin + end) >> 1;
        uint64_t a_val = A[mid];
        uint64_t b_val = B[diag - 1 - mid];
        if (a_val < b_val) {
            begin = mid + 1;
        } else {
            end = mid;
        }
    }
    return begin;
}

// 2. Kernel de Fusão Paralela
__global__ void merge_path_kernel(const uint64_t* A, int lenA, const uint64_t* B, int lenB, uint64_t* C) {
    int total_elements = lenA + lenB;
    int elements_per_thread = 8;
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

    int diag_start = min(thread_id * elements_per_thread, total_elements);
    int diag_end = min(diag_start + elements_per_thread, total_elements);

    if (diag_start >= total_elements) return;

    int start_A = merge_path_search(A, lenA, B, lenB, diag_start);
    int start_B = diag_start - start_A;
    int end_A = merge_path_search(A, lenA, B, lenB, diag_end);
    int end_B = diag_end - end_A;

    int idx_A = start_A;
    int idx_B = start_B;
    int idx_C = diag_start;

    while (idx_A < end_A && idx_B < end_B) {
        if (A[idx_A] <= B[idx_B]) { C[idx_C++] = A[idx_A++]; }
        else { C[idx_C++] = B[idx_B++]; }
    }
    while (idx_A < end_A) { C[idx_C++] = A[idx_A++]; }
    while (idx_B < end_B) { C[idx_C++] = B[idx_B++]; }
}

// 3. FFI Wrapper COM TRANSFERÊNCIA DE MEMÓRIA (O Segredo do Zero-Copy Real)
extern "C" {
    int launch_merge_path_kernel(const uint64_t* h_A, const uint64_t* h_B, uint64_t* h_C, size_t len_A, size_t len_B) {
        uint64_t *d_A = NULL, *d_B = NULL, *d_C = NULL;
        size_t size_A = len_A * sizeof(uint64_t);
        size_t size_B = len_B * sizeof(uint64_t);
        size_t size_C = (len_A + len_B) * sizeof(uint64_t);
        int erro = 0;
        int total_elements = 0;
        int blocks = 0;
        const int elements_per_thread = 8;
        const int threads_per_block = 256;

        // O kernel indexa com `int`: mais que isso ele não endereça.
        if (len_A + len_B > (size_t)0x7fffffff) return MM_ERRO_TAMANHO;
        if (len_A + len_B == 0) return 0;

        // 1. Aloca memória na VRAM da NVIDIA (1 byte a mais evita pedido vazio)
        MM_TENTA(cudaMalloc((void**)&d_A, size_A + 1));
        MM_TENTA(cudaMalloc((void**)&d_B, size_B + 1));
        MM_TENTA(cudaMalloc((void**)&d_C, size_C));

        // 2. Copia os dados da Placa-Mãe (RAM) para a Placa de Vídeo (VRAM)
        if (size_A) MM_TENTA(cudaMemcpy(d_A, h_A, size_A, cudaMemcpyHostToDevice));
        if (size_B) MM_TENTA(cudaMemcpy(d_B, h_B, size_B, cudaMemcpyHostToDevice));

        // 3. Executa o Kernel C++
        total_elements = (int)(len_A + len_B);
        blocks = (int)(((size_t)total_elements + (size_t)(elements_per_thread * threads_per_block) - 1) / (size_t)(elements_per_thread * threads_per_block));
        if (blocks == 0) blocks = 1;

        merge_path_kernel<<<blocks, threads_per_block>>>(d_A, (int)len_A, d_B, (int)len_B, d_C);
        MM_TENTA(cudaGetLastError());
        MM_TENTA(cudaDeviceSynchronize());

        // 4. Copia o resultado final de volta para o Rust
        MM_TENTA(cudaMemcpy(h_C, d_C, size_C, cudaMemcpyDeviceToHost));

    fim:
        // 5. Limpa a GPU (cudaFree aceita NULL)
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
        if (erro) cudaGetLastError();
        return erro;
    }
}
// =========================================================
// MÓDULO 6: FILTRO BLOOM PROBABILÍSTICO (CUDA)
// =========================================================
// ATENCAO: o hash, as sementes e o arranjo dos bits abaixo tem um gemeo na CPU
// (src/bloom.rs). Mudar um sem mudar o outro faz o filtro gravado por um nao
// ser lido pelo outro.

// Função criptográfica para GPU (MurmurHash-style)
__device__ uint64_t gpu_hash_kmer(uint64_t k, uint32_t seed) {
    k ^= seed;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

// Kernel que roda em milhares de núcleos simultaneamente
__global__ void bloom_insert_kernel(const uint64_t* kmers, size_t num_kmers, uint32_t* bloom, size_t bloom_bits, int num_hashes) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_kmers) return;

    uint64_t kmer = kmers[tid];

    // Calcula múltiplos Hashes para o mesmo K-mer e "acende" os bits
    for (int i = 0; i < num_hashes; i++) {
        uint64_t hash_val = gpu_hash_kmer(kmer, (uint32_t)i * 0x9E3779B9u);
        size_t bit_idx = hash_val % bloom_bits;

        size_t word_idx = bit_idx / 32;
        uint32_t bit_mask = 1u << (bit_idx % 32);

        // atomicOr previne que duas threads sobrescrevam a memória ao mesmo tempo!
        atomicOr(&bloom[word_idx], bit_mask);
    }
}

extern "C" {
    // A Ponte C++ para o Rust
    int launch_bloom_filter_build(const uint64_t* h_kmers, size_t num_kmers, uint8_t* h_bloom, size_t bloom_bytes, int num_hashes) {
        uint64_t* d_kmers = NULL;
        uint32_t* d_bloom = NULL;
        size_t bloom_words = bloom_bytes / 4;
        size_t bloom_bits = bloom_words * 32;
        int erro = 0;
        const int threads = 256;
        size_t blocks = (num_kmers + threads - 1) / threads;

        // Sem bits nao ha onde gravar; sem k-mers o filtro e todo zero.
        if (bloom_bits == 0 || blocks > (size_t)0x7fffffff) return MM_ERRO_TAMANHO;
        if (num_kmers == 0) {
            memset(h_bloom, 0, bloom_bytes);
            return 0;
        }

        // 1. Aloca memória na VRAM
        MM_TENTA(cudaMalloc((void**)&d_kmers, num_kmers * sizeof(uint64_t)));
        MM_TENTA(cudaMalloc((void**)&d_bloom, bloom_bytes));

        // 2. Transfere os genomas e zera o Filtro
        MM_TENTA(cudaMemcpy(d_kmers, h_kmers, num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice));
        MM_TENTA(cudaMemset(d_bloom, 0, bloom_bytes));

        // 3. Força Bruta: Lança milhares de threads!
        bloom_insert_kernel<<<(unsigned int)blocks, threads>>>(d_kmers, num_kmers, d_bloom, bloom_bits, num_hashes);
        MM_TENTA(cudaGetLastError());
        MM_TENTA(cudaDeviceSynchronize());

        // 4. Devolve o Filtro compactado para a Memória RAM
        MM_TENTA(cudaMemcpy(h_bloom, d_bloom, bloom_bytes, cudaMemcpyDeviceToHost));

    fim:
        cudaFree(d_kmers);
        cudaFree(d_bloom);
        if (erro) cudaGetLastError();
        return erro;
    }
}
// =========================================================
// MÓDULO 4: QUERY DE VARIANTES (SNPs) NO FILTRO BLOOM
// =========================================================

__global__ void bloom_query_kernel(const uint64_t* query_kmers, size_t num_queries, const uint32_t* bloom, size_t bloom_bits, int num_hashes, uint8_t* out_mask) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;

    uint64_t kmer = query_kmers[tid];
    uint8_t exists = 1; // Assume que o K-mer é saudável (existe na referência)

    for (int i = 0; i < num_hashes; i++) {
        uint64_t hash_val = gpu_hash_kmer(kmer, (uint32_t)i * 0x9E3779B9u);
        size_t bit_idx = hash_val % bloom_bits;
        size_t word_idx = bit_idx / 32;
        uint32_t bit_mask = 1u << (bit_idx % 32);

        // Se encontrar 1 único bit desligado, é impossível que este DNA esteja na referência! É uma mutação!
        if ((bloom[word_idx] & bit_mask) == 0) {
            exists = 0;
            break;
        }
    }
    out_mask[tid] = exists;
}

extern "C" {
    int launch_bloom_filter_query(const uint64_t* h_query_kmers, size_t num_queries, const uint8_t* h_bloom, size_t bloom_bytes, int num_hashes, uint8_t* h_out_mask) {
        uint64_t* d_query_kmers = NULL;
        uint32_t* d_bloom = NULL;
        uint8_t* d_out_mask = NULL;
        size_t bloom_words = bloom_bytes / 4;
        size_t bloom_bits = bloom_words * 32;
        int erro = 0;
        const int threads = 256;
        size_t blocks = (num_queries + threads - 1) / threads;

        if (bloom_bits == 0 || blocks > (size_t)0x7fffffff) return MM_ERRO_TAMANHO;
        if (num_queries == 0) return 0;

        // 1. Aloca memória na VRAM
        MM_TENTA(cudaMalloc((void**)&d_query_kmers, num_queries * sizeof(uint64_t)));
        MM_TENTA(cudaMalloc((void**)&d_bloom, bloom_bytes));
        MM_TENTA(cudaMalloc((void**)&d_out_mask, num_queries * sizeof(uint8_t)));

        // 2. Transfere os dados (Amostra e o Filtro já construído)
        MM_TENTA(cudaMemcpy(d_query_kmers, h_query_kmers, num_queries * sizeof(uint64_t), cudaMemcpyHostToDevice));
        MM_TENTA(cudaMemcpy(d_bloom, h_bloom, bloom_bytes, cudaMemcpyHostToDevice));

        // 3. Força Bruta: Lança milhares de threads para inspecionar os K-mers!
        bloom_query_kernel<<<(unsigned int)blocks, threads>>>(d_query_kmers, num_queries, d_bloom, bloom_bits, num_hashes, d_out_mask);
        MM_TENTA(cudaGetLastError());
        MM_TENTA(cudaDeviceSynchronize());

        // 4. Devolve a máscara de resultados (0 = Mutação, 1 = Saudável)
        MM_TENTA(cudaMemcpy(h_out_mask, d_out_mask, num_queries * sizeof(uint8_t), cudaMemcpyDeviceToHost));

    fim:
        cudaFree(d_query_kmers);
        cudaFree(d_bloom);
        cudaFree(d_out_mask);
        if (erro) cudaGetLastError();
        return erro;
    }
}
