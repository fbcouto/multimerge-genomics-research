#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

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
    void launch_merge_path_kernel(const uint64_t* h_A, const uint64_t* h_B, uint64_t* h_C, size_t len_A, size_t len_B) {
        uint64_t *d_A, *d_B, *d_C;
        size_t size_A = len_A * sizeof(uint64_t);
        size_t size_B = len_B * sizeof(uint64_t);
        size_t size_C = (len_A + len_B) * sizeof(uint64_t);

        // 1. Aloca memória na VRAM da NVIDIA
        cudaMalloc((void**)&d_A, size_A);
        cudaMalloc((void**)&d_B, size_B);
        cudaMalloc((void**)&d_C, size_C);

        // 2. Copia os dados da Placa-Mãe (RAM) para a Placa de Vídeo (VRAM)
        cudaMemcpy(d_A, h_A, size_A, cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, h_B, size_B, cudaMemcpyHostToDevice);

        // 3. Executa o Kernel C++
        int total_elements = len_A + len_B;
        int elements_per_thread = 8;
        int threads_per_block = 256;
        int blocks = (total_elements + (elements_per_thread * threads_per_block) - 1) / (elements_per_thread * threads_per_block);
        if (blocks == 0) blocks = 1;

        merge_path_kernel<<<blocks, threads_per_block>>>(d_A, (int)len_A, d_B, (int)len_B, d_C);
        cudaDeviceSynchronize();

        // 4. Copia o resultado final de volta para o Rust
        cudaMemcpy(h_C, d_C, size_C, cudaMemcpyDeviceToHost);

        // 5. Limpa a GPU
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }
}
// =========================================================
// MÓDULO 6: FILTRO BLOOM PROBABILÍSTICO (CUDA)
// =========================================================

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
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_kmers) return;

    uint64_t kmer = kmers[tid];
    
    // Calcula múltiplos Hashes para o mesmo K-mer e "acende" os bits
    for (int i = 0; i < num_hashes; i++) {
        uint64_t hash_val = gpu_hash_kmer(kmer, i * 0x9E3779B9);
        size_t bit_idx = hash_val % bloom_bits;
        
        size_t word_idx = bit_idx / 32;
        uint32_t bit_mask = 1 << (bit_idx % 32);
        
        // atomicOr previne que duas threads sobrescrevam a memória ao mesmo tempo!
        atomicOr(&bloom[word_idx], bit_mask);
    }
}

extern "C" {
    // A Ponte C++ para o Rust
    void launch_bloom_filter_build(const uint64_t* h_kmers, size_t num_kmers, uint8_t* h_bloom, size_t bloom_bytes, int num_hashes) {
        uint64_t* d_kmers;
        uint32_t* d_bloom;
        size_t bloom_words = bloom_bytes / 4;
        size_t bloom_bits = bloom_words * 32;

        // 1. Aloca memória na VRAM
        cudaMalloc((void**)&d_kmers, num_kmers * sizeof(uint64_t));
        cudaMalloc((void**)&d_bloom, bloom_bytes);

        // 2. Transfere os genomas e zera o Filtro
        cudaMemcpy(d_kmers, h_kmers, num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemset(d_bloom, 0, bloom_bytes);

        // 3. Força Bruta: Lança milhares de threads!
        int threads = 256;
        int blocks = (num_kmers + threads - 1) / threads;
        bloom_insert_kernel<<<blocks, threads>>>(d_kmers, num_kmers, d_bloom, bloom_bits, num_hashes);
        cudaDeviceSynchronize();

        // 4. Devolve o Filtro compactado para a Memória RAM
        cudaMemcpy(h_bloom, d_bloom, bloom_bytes, cudaMemcpyDeviceToHost);

        cudaFree(d_kmers);
        cudaFree(d_bloom);
    }
}
// =========================================================
// MÓDULO 4: QUERY DE VARIANTES (SNPs) NO FILTRO BLOOM
// =========================================================

__global__ void bloom_query_kernel(const uint64_t* query_kmers, size_t num_queries, const uint32_t* bloom, size_t bloom_bits, int num_hashes, uint8_t* out_mask) {
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;

    uint64_t kmer = query_kmers[tid];
    uint8_t exists = 1; // Assume que o K-mer é saudável (existe na referência)
    
    for (int i = 0; i < num_hashes; i++) {
        uint64_t hash_val = gpu_hash_kmer(kmer, i * 0x9E3779B9);
        size_t bit_idx = hash_val % bloom_bits;
        size_t word_idx = bit_idx / 32;
        uint32_t bit_mask = 1 << (bit_idx % 32);
        
        // Se encontrar 1 único bit desligado, é impossível que este DNA esteja na referência! É uma mutação!
        if ((bloom[word_idx] & bit_mask) == 0) {
            exists = 0; 
            break;
        }
    }
    out_mask[tid] = exists;
}

extern "C" {
    void launch_bloom_filter_query(const uint64_t* h_query_kmers, size_t num_queries, const uint8_t* h_bloom, size_t bloom_bytes, int num_hashes, uint8_t* h_out_mask) {
        uint64_t* d_query_kmers;
        uint32_t* d_bloom;
        uint8_t* d_out_mask;
        size_t bloom_words = bloom_bytes / 4;
        size_t bloom_bits = bloom_words * 32;

        // 1. Aloca memória na VRAM
        cudaMalloc((void**)&d_query_kmers, num_queries * sizeof(uint64_t));
        cudaMalloc((void**)&d_bloom, bloom_bytes);
        cudaMalloc((void**)&d_out_mask, num_queries * sizeof(uint8_t));

        // 2. Transfere os dados (Amostra e o Filtro já construído)
        cudaMemcpy(d_query_kmers, h_query_kmers, num_queries * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_bloom, h_bloom, bloom_bytes, cudaMemcpyHostToDevice);

        // 3. Força Bruta: Lança milhares de threads para inspecionar os K-mers!
        int threads = 256;
        int blocks = (num_queries + threads - 1) / threads;
        bloom_query_kernel<<<blocks, threads>>>(d_query_kmers, num_queries, d_bloom, bloom_bits, num_hashes, d_out_mask);
        cudaDeviceSynchronize();

        // 4. Devolve a máscara de resultados (0 = Mutação, 1 = Saudável)
        cudaMemcpy(h_out_mask, d_out_mask, num_queries * sizeof(uint8_t), cudaMemcpyDeviceToHost);

        cudaFree(d_query_kmers);
        cudaFree(d_bloom);
        cudaFree(d_out_mask);
    }
}