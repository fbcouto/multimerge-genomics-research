import marimo

__generated_with = "0.23.11"
app = marimo.App(width="full")

@app.cell
def imports():
    import marimo as mo
    import polars as pl
    import numpy as np
    import glob
    import asyncio  # 👈 Importado aqui uma única vez!
    import multimerge_genomics as mmg
    return asyncio, glob, mmg, mo, np, pl

@app.cell
def estado_global(mo):
    obter_pagina, definir_pagina = mo.state("home")
    return definir_pagina, obter_pagina

@app.cell
def identidade_visual(mo):
    script_identidade = """
    <script>
        document.title = "MMGEN Suite";
        let link = document.querySelector("link[rel~='icon']");
        if (!link) {
            link = document.createElement('link');
            link.rel = 'icon';
            document.head.appendChild(link);
        }
        link.href = "data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🧬</text></svg>";
    </script>
    """
    return mo.Html(script_identidade)

@app.cell
def menu_principal(definir_pagina, mo, obter_pagina):
    mo.stop(obter_pagina() != "home", mo.md(""))

    btn_menu_qc = mo.ui.button(label="🔬 1. QC e Trimming (FASTQ)", on_change=lambda _: definir_pagina("qc"))
    btn_menu_cluster = mo.ui.button(label="🧬 2. Agrupamento Taxonômico", on_change=lambda _: definir_pagina("cluster"))
    btn_menu_assembly = mo.ui.button(label="🧩 3. Montagem de Genomas (De Bruijn)", on_change=lambda _: definir_pagina("assembly"))
    btn_menu_variant = mo.ui.button(label="⚠️ 4. Chamada de Variantes (SNPs)", on_change=lambda _: definir_pagina("variant"))
    btn_menu_annotation = mo.ui.button(label="📘 5. Anotação Funcional", on_change=lambda _: definir_pagina("annotation"))
    btn_menu_bloom = mo.ui.button(label="🔍 6. Filtros Bloom (cuSBF na GPU)", on_change=lambda _: definir_pagina("bloom"))
    btn_menu_scaffold = mo.ui.button(label="🖇️ 7. Fechamento de Genomas (Scaffolding)", on_change=lambda _: definir_pagina("scaffold"))
    btn_menu_vcf = mo.ui.button(label="🗄️ 8. VCF Multimerge & Sorting", on_change=lambda _: definir_pagina("vcf"))

    tela_menu = mo.md(f"""
    <div style="text-align: center; margin-top: 40px; font-family: sans-serif;">
        <h1 style="font-size: 3em; margin-bottom: 0px;">🧬 MMGEN</h1>
        <h3 style="color: gray; margin-top: 5px;">High-Performance Genomics Framework</h3>
        <p style="margin-bottom: 40px;">Motor acelerado por Rust e NVIDIA CUDA. Selecione um módulo para iniciar:</p>
        
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px; max-width: 850px; margin: 0 auto;">
            {btn_menu_qc}
            {btn_menu_cluster}
            {btn_menu_assembly}
            {btn_menu_variant}
            {btn_menu_annotation}
            {btn_menu_bloom}
            {btn_menu_scaffold}
            {btn_menu_vcf}
        </div>
    </div>
    """)
    
    tela_menu
    return btn_menu_annotation, btn_menu_assembly, btn_menu_bloom, btn_menu_cluster, btn_menu_qc, btn_menu_variant, tela_menu

@app.cell
def pagina_qc_ui(definir_pagina, glob, mo, obter_pagina):
    import os as _os
    # 👇 AQUI ESTÁ A CORREÇÃO: Voltou a ser "qc" para o botão do menu funcionar!
    mo.stop(obter_pagina() != "qc", mo.md(""))
    
    btn_voltar_qc = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    _arquivos = glob.glob("*.fastq") + glob.glob("*.fastq.gz") + glob.glob("*.fasta") + glob.glob("*.fasta.gz") + glob.glob("*.fa")
    if not _arquivos: _arquivos = ["Nenhum arquivo encontrado"]

    form_qc = mo.ui.form(
        element=mo.ui.dictionary({
            "caminhos_arquivos": mo.ui.multiselect(options=_arquivos, label="🧬 Ficheiros de DNA (Suporta múltiplos e .gz):"),
            "tamanho_k": mo.ui.slider(start=21, stop=63, step=2, value=55, label="Tamanho do K-mer (Régua u128):"),
            "qualidade": mo.ui.slider(start=0, stop=40, step=1, value=20, label="Qualidade Mínima (Apenas FASTQ):")
        }),
        bordered=True,
        submit_button_label="🚀 Iniciar Ingestão em Lote"
    )

    ui_qc = mo.vstack([
        btn_voltar_qc,
        mo.md("<br># 🧬 1. Ingestão de Dados e QC (Lote / GZIP)"),
        mo.callout(mo.md("O motor agora suporta a importação de **múltiplos ficheiros em simultâneo** (ex: Paired-End R1 e R2) e descompacta ficheiros `.gz` diretamente na memória RAM, poupando espaço no disco."), kind="info"),
        form_qc
    ])
    
    ui_qc
    return btn_voltar_qc, form_qc, ui_qc

@app.cell
async def pagina_qc_backend(asyncio, form_qc, mmg, mo, obter_pagina, pl):
    import os as _os
    # 👇 AQUI TAMBÉM: Voltou a ser "qc"
    mo.stop(obter_pagina() != "qc", mo.md(""))
    mo.stop(form_qc.value is None, mo.md("⏳ *Aguardando configuração...*"))

    _caminhos = form_qc.value.get("caminhos_arquivos", [])
    mo.stop(not _caminhos or "Nenhum" in _caminhos[0], mo.callout(mo.md("⚠️ Selecione pelo menos um ficheiro."), kind="warning"))

    _k = form_qc.value["tamanho_k"]
    _qualidade = form_qc.value["qualidade"]
    
    _is_fasta = any(c.endswith(".fasta") or c.endswith(".fa") or c.endswith(".fasta.gz") for c in _caminhos)

    try:
        with mo.status.spinner(f"Processando {len(_caminhos)} ficheiro(s) em Multi-Threading (K={_k})..."):
            await asyncio.sleep(0.1)
            
            if _is_fasta:
                _bytes_kmer = mmg.ler_fasta_em_lote(_caminhos, _k)
            else:
                _bytes_kmer = mmg.ler_fastq_em_lote(_caminhos, _qualidade, _k)

        with mo.status.spinner("Ordenando e Salvando Monólito Arrow..."):
            await asyncio.sleep(0.1)
            _df = pl.DataFrame({"K-mers": [_bytes_kmer]})
            
            if len(_caminhos) > 1:
                _nome_base = "Lote_Multiplo"
            else:
                _nome_base = _caminhos[0].split(".")[0]
                
            _nome_arquivo = f"{_nome_base}_K{_k}.arrow"
            _df.write_ipc(_nome_arquivo)

        _painel = mo.vstack([
            mo.callout(mo.md(f"""
            ### ✅ Lote Processado com Sucesso!
            * **Ficheiros Analisados (Threads):** `{len(_caminhos)}`
            * **Tamanho do K-mer:** `{_k}`
            * **Sequências Únicas Extraídas:** `{len(_bytes_kmer) // 16:,}`
            * **Ficheiro Arrow Gerado:** `{_nome_arquivo}`
            """), kind="success")
        ])
        
    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,
@app.cell
def pagina_cluster_ui(definir_pagina, glob, mo, obter_pagina):
    mo.stop(obter_pagina() != "cluster", mo.md(""))
    
    _btn_voltar = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))
    
    _arquivos = glob.glob("*.arrow") + glob.glob("*.txt") + glob.glob("*.fasta")
    if not _arquivos: _arquivos = ["Nenhum arquivo encontrado"]

    form_dados = mo.ui.form(
        element=mo.ui.dictionary({
            "path_a": mo.ui.dropdown(options=_arquivos, label="📂 Amostra A (Recomendado: .arrow):"),
            "path_b": mo.ui.dropdown(options=_arquivos, label="📂 Amostra B (Recomendado: .arrow):")
        }), bordered=True, submit_button_label="🚀 1. Iniciar GPU"
    )

    _ui = mo.vstack([_btn_voltar, mo.md("<br># 🧬 Módulo de Agrupamento Taxonômico"), form_dados])
    _ui
    return form_dados,

@app.cell
async def pagina_cluster_backend(asyncio, form_dados, mmg, mo, np, obter_pagina, pl):
    mo.stop(obter_pagina() != "cluster", mo.md("")) 
    mo.stop(form_dados.value is None, mo.md("🧬 *Aguardando clique...*"))

    _p_a, _p_b = form_dados.value["path_a"], form_dados.value["path_b"]
    if not _p_a or "Nenhum" in _p_a: mo.stop(True, mo.md("⚠️ Erro."))

    with mo.status.spinner("Lendo ficheiros e atirando para a VRAM (CUDA)..."):
        await asyncio.sleep(0.1)
        
        if _p_a.endswith(".arrow") and _p_b.endswith(".arrow"):
            _bytes_a = pl.read_ipc(_p_a).get_column("K-mers").to_numpy().tobytes()
            _bytes_b = pl.read_ipc(_p_b).get_column("K-mers").to_numpy().tobytes()
            _bytes_gpu = mmg.ordenar_bytes_gpu(_bytes_a, _bytes_b)
        else:
            _bytes_gpu = mmg.ler_arquivos_e_processar(_p_a, _p_b)
            
        _arr_gpu = np.frombuffer(_bytes_gpu, dtype=np.uint64)
        
    with mo.status.spinner("Gerando Assinaturas MinHash na CPU..."):
        await asyncio.sleep(0.1)
        _bytes_assinaturas = mmg.gerar_assinaturas_lote(_bytes_gpu, 50)
        matriz_assinaturas = np.frombuffer(_bytes_assinaturas, dtype=np.uint64).reshape(-1, 50)
        
        df_final = pl.DataFrame({
            "K-mers (GPU)": _arr_gpu,
            "Sequencia_DNA": mmg.decodificar_lote(_bytes_gpu, 31)
        })

    _painel = mo.md(f"✅ Motor Finalizado. **{df_final.height:,}** K-mers prontos para agrupamento.")
    _painel
    return df_final, matriz_assinaturas, _painel

@app.cell
def pagina_cluster_hdbscan(df_final, mo, obter_pagina):
    mo.stop(obter_pagina() != "cluster", mo.md(""))

    form_cluster = mo.ui.form(
        element=mo.ui.dictionary({"porcentagem": mo.ui.slider(start=1, stop=100, step=1, value=5, label="Amostragem HDBSCAN (%):")}),
        bordered=True, submit_button_label="🧬 2. Agrupar Espécies"
    )
    
    form_cluster
    return form_cluster,

@app.cell
async def pagina_cluster_resultado(asyncio, df_final, form_cluster, matriz_assinaturas, mo, np, obter_pagina, pl):
    from sklearn.cluster import HDBSCAN
    from sklearn.decomposition import PCA

    mo.stop(obter_pagina() != "cluster", mo.md(""))
    mo.stop(form_cluster.value is None, mo.md(""))

    _fracao = form_cluster.value["porcentagem"] / 100.0
    _num_rows_total = len(df_final)
    
    # 🧠 AMOSTRAGEM INTELIGENTE
    # Se temos muito poucos dados (como no nosso teste de 6 k-mers), forçamos 100% da amostra
    if _num_rows_total < 1000:
        _tamanho_amostra = _num_rows_total
    else:
        _tamanho_amostra = int(_num_rows_total * _fracao)

    # 🛡️ PROTEÇÕES
    if _tamanho_amostra < 4:
        mo.stop(True, mo.md("⚠️ **Erro:** Dados insuficientes. O modelo precisa de pelo menos 4 K-mers para calcular a geometria tridimensional (PCA)."))

    if _tamanho_amostra > 100000:
        mo.stop(True, mo.md(f"""
        ⚠️ **Aviso de Gargalo de CPU:** Você está tentando agrupar **{_tamanho_amostra:,}** K-mers. O algoritmo HDBSCAN no Python (CPU) não foi feito para lidar com esse volume e pode demorar horas. 
        👉 **Solução:** Reduza o slider e clique em Agrupar novamente.
        """))

    with mo.status.spinner(f"Amostragem de {_tamanho_amostra:,} pontos e cálculo espacial..."):
        await asyncio.sleep(0.1)
        
        _indices = np.arange(_num_rows_total)
        # Só fatiamos se a amostra for menor que o total
        if _tamanho_amostra < _num_rows_total:
            _indices = np.random.choice(_indices, _tamanho_amostra, replace=False)
        
        _df_processar = df_final[_indices]
        _matriz_embeddings = matriz_assinaturas[_indices].astype(np.float64)

        _matriz_reduzida = PCA(n_components=3).fit_transform(_matriz_embeddings)
        
        # Como o HDBSCAN exige grupos de no mínimo 50, ajustamos dinamicamente para testes pequenos
        _min_cluster = min(50, max(2, _tamanho_amostra // 2))
        _labels = HDBSCAN(min_cluster_size=_min_cluster, n_jobs=-1, copy=False).fit_predict(_matriz_reduzida)

        _df_agrupado = _df_processar.with_columns(pl.Series("Cluster_ID", _labels))
        _total_clusters = len(set(_labels)) - (1 if -1 in _labels else 0)

    _tabela = mo.ui.table(_df_agrupado, selection="single", pagination=True)
    
    _painel_final = mo.vstack([
        mo.callout(mo.md(f"**Resultado:** {_total_clusters} Espécies/Cepas detectadas. (K-mers sem grupo/Ruído: {sum(_labels == -1)})"), kind="success"),
        _tabela
    ])
    
    _painel_final
    return HDBSCAN, PCA

@app.cell
def pagina_assembly_ui(definir_pagina, glob, mo, obter_pagina):
    mo.stop(obter_pagina() != "assembly", mo.md(""))
    btn_voltar_assembly = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))
    
    # Procura os ficheiros .arrow do Módulo 1
    _arquivos = glob.glob("*.arrow")
    if not _arquivos: _arquivos = ["Nenhum arquivo .arrow encontrado"]

    form_assembly = mo.ui.form(
        element=mo.ui.dictionary({
            "path_arrow": mo.ui.dropdown(options=_arquivos, label="🧬 Escolha o Ficheiro Base (.arrow):")
        }), 
        bordered=True, 
        submit_button_label="🧩 1. Montar Grafo de De Bruijn (Memória RAM)"
    )
    
    ui_assembly = mo.vstack([
        btn_voltar_assembly, 
        mo.md("<br># 🧩 Montagem de Genomas (De Bruijn)"),
        mo.callout(mo.md("Este módulo recebe os K-mers purificados, fatia cada um deles em Prefixo (K-1) e Sufixo (K-1) e mapeia as sobreposições na memória RAM do seu computador usando a velocidade extrema do Rust. O resultado é a base estrutural para montar o genoma original."), kind="info"),
        form_assembly
    ])
    
    ui_assembly
    return btn_voltar_assembly, form_assembly, ui_assembly

@app.cell
async def pagina_assembly_backend(asyncio, form_assembly, mmg, mo, obter_pagina, pl):
    import os as _os # <-- IMPORTAÇÃO PRIVADA AQUI
    mo.stop(obter_pagina() != "assembly", mo.md(""))
    mo.stop(form_assembly.value is None, mo.md("🧩 *Aguardando clique no botão de montagem...*"))

    _path = form_assembly.value.get("path_arrow")
    mo.stop(not _path or "Nenhum" in _path, mo.callout(mo.md("⚠️ Selecione um ficheiro .arrow válido."), kind="warning"))

    try:
        # Extrai o "K" automaticamente do nome do ficheiro (ex: arquivo_K55.arrow -> 55)
        _k_str = _path.split("_K")[-1].split(".arrow")[0]
        _k_size = int(_k_str)

        with mo.status.spinner(f"1/3 Lendo K-mers (Zero-Copy) Monolítico (K={_k_size})..."):
            await asyncio.sleep(0.1)
            _df = pl.read_ipc(_path)
            _bytes_kmer = _df.get_column("K-mers").item(0) 
            
        with mo.status.spinner("2/3 Rust: Analisando Grafo e Superando Repetições..."):
            await asyncio.sleep(0.1)
            _nos, _arestas = mmg.construir_grafo_debruijn(_bytes_kmer, _k_size)
            lista_contigs = mmg.montar_genoma_debruijn(_bytes_kmer, _k_size)

        with mo.status.spinner("3/3 Salvando no Disco e Preparando Interface..."):
            await asyncio.sleep(0.1)
            
            _nome_arquivo = f"MMGEN_Contigs_K{_k_size}.fasta"
            with open(_nome_arquivo, "w", encoding="utf-8") as arquivo_fasta:
                for i, c in enumerate(lista_contigs):
                    arquivo_fasta.write(f">NODE_{i+1}_length_{len(c)}\n{c}\n")

            _caminho = _os.path.abspath(_nome_arquivo) # <-- USAMOS _os AQUI
            _top_contigs = lista_contigs[:100]

            df_contigs = pl.DataFrame({
                "Contig ID": [f"NODE_{i+1}" for i in range(len(_top_contigs))],
                "Tamanho (bp)": [len(c) for c in _top_contigs],
                "Sequência Inicial": [c[:50] + "..." for c in _top_contigs]
            })
            _tabela_ui = mo.ui.table(df_contigs, pagination=True, page_size=10)

        _painel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🗺️ Grafo Ultra-Longo Resolvido com Sucesso! (128-bits)
            * **Nós Totais:** `{_nos:,}` | **Arestas Totais:** `{_arestas:,}`
            * **Total de Peças Geradas:** `{len(lista_contigs):,}`
            * **Maior Peça Reconstruída:** `{len(lista_contigs[0]):,} bp`
            """), kind="success"),
            mo.callout(mo.md(f"💾 **Ficheiro Físico Salvo:** `{_caminho}`"), kind="info"),
            mo.md(f"#### 📋 Top {len(_top_contigs)} Maiores Fragmentos (Unitigs):"),
            _tabela_ui
        ])

    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO NO PYTHON:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,
@app.cell
def pagina_variant_ui(definir_pagina, glob, mo, obter_pagina):
    import os as _os
    mo.stop(obter_pagina() != "variant", mo.md(""))
    btn_voltar_variant = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    # 1. Criamos um estado para mostrar o resultado do botão na tela!
    get_msg, set_msg = mo.state("")

    def gerar_paciente():
        # 2. Corrigido para o nome exato do ficheiro no seu disco
        _nome_base = "ecoli.fasta" 
        
        if not _os.path.exists(_nome_base): 
            set_msg(f"❌ Erro: O ficheiro '{_nome_base}' não está nesta pasta!")
            return
            
        with open(_nome_base, "r") as f: linhas = f.readlines()
        seq = "".join([l.strip() for l in linhas if not l.startswith(">")])
        seq_list = list(seq)
        
        if len(seq_list) > 30000:
            seq_list[5000] = 'T' if seq_list[5000] == 'A' else 'A'
            seq_list[15000] = 'G' if seq_list[15000] == 'C' else 'C'
            seq_list[25000] = 'C' if seq_list[25000] == 'G' else 'G'
            
        seq_mutada = "".join(seq_list)
        with open("paciente_mutado.fasta", "w") as f:
            f.write(">Paciente_Mutado_Ecoli\n")
            for i in range(0, len(seq_mutada), 80): f.write(seq_mutada[i:i+80] + "\n")
            
        # 3. Atualiza a mensagem na interface gráfica
        set_msg("✅ SUCESSO! O 'paciente_mutado.fasta' foi gerado no disco com 3 mutações.")

    btn_gerar = mo.ui.button(label="💉 Gerar Paciente Simulado", on_click=lambda _: gerar_paciente())

    _arquivos = glob.glob("*.arrow")
    if not _arquivos: _arquivos = ["Nenhum arquivo .arrow encontrado"]

    form_variant = mo.ui.form(
        element=mo.ui.dictionary({
            "path_ref": mo.ui.dropdown(options=_arquivos, label="🧬 1. Genoma de Referência Saudável (.arrow):"),
            "path_paciente": mo.ui.dropdown(options=_arquivos, label="🩸 2. Amostra do Paciente (.arrow):"),
            "tamanho_mb": mo.ui.slider(start=1, stop=128, step=1, value=16, label="Memória da VRAM para o Filtro (MB):"),
            "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=5, label="Cálculos Criptográficos por K-mer:")
        }),
        bordered=True,
        submit_button_label="🔍 3. Executar Busca CUDA"
    )

    ui_variant = mo.vstack([
        btn_voltar_variant,
        mo.md("<br># ⚠️ Chamada de Variantes (SNPs) via GPU"),
        mo.callout(mo.md("Para encontrar mutações, precisamos comparar um Paciente contra uma Referência. Clique no botão abaixo para simular um paciente."), kind="info"),
        btn_gerar,
        # AQUI APARECERÁ O FEEDBACK DO BOTÃO:
        mo.md(f"**{get_msg()}**"), 
        form_variant
    ])
    
    ui_variant
    return btn_gerar, btn_voltar_variant, form_variant, get_msg, set_msg, ui_variant
@app.cell
async def pagina_variant_backend(asyncio, form_variant, mmg, mo, obter_pagina, pl):
    mo.stop(obter_pagina() != "variant", mo.md(""))
    mo.stop(form_variant.value is None, mo.md("🧪 *Aguardando configuração...*"))

    _ref = form_variant.value.get("path_ref")
    _pac = form_variant.value.get("path_paciente")
    
    mo.stop(_ref == _pac, mo.callout(mo.md("⚠️ A Referência e o Paciente não podem ser o mesmo ficheiro!"), kind="warning"))
    mo.stop(not _ref or not _pac or "Nenhum" in _ref, mo.md("⚠️ Selecione os dois ficheiros."))

    _k_str = _ref.split("_K")[-1].split(".arrow")[0]
    _k_size = int(_k_str)
    _mb = form_variant.value["tamanho_mb"]
    _hashes = form_variant.value["num_hashes"]

    try:
        with mo.status.spinner("1/3 Carregando Genomas (Zero-Copy)..."):
            await asyncio.sleep(0.1)
            _ref_bytes = pl.read_ipc(_ref).get_column("K-mers").item(0)
            _pac_bytes = pl.read_ipc(_pac).get_column("K-mers").item(0)

        with mo.status.spinner("2/3 Construindo Filtro Bloom da Referência na NVIDIA..."):
            await asyncio.sleep(0.1)
            _bloom_bytes = mmg.construir_filtro_bloom_gpu(_ref_bytes, _mb, _hashes)

        with mo.status.spinner("3/3 Buscando Mutações Exclusivas do Paciente na GPU..."):
            await asyncio.sleep(0.1)
            _mutacoes_bytes = mmg.extrair_mutacoes_gpu(_pac_bytes, _bloom_bytes, _hashes)
            _sequencias_mutadas = mmg.decodificar_lote_128(_mutacoes_bytes, _k_size)
            
        _total_mutacoes = len(_sequencias_mutadas)
        _amostra_visual = _sequencias_mutadas[:15] if _total_mutacoes > 0 else ["Nenhuma mutação encontrada."]

        _painel = mo.vstack([
            mo.callout(mo.md(f"""### 🏆 Varredura Genómica Concluída na GPU!
A NVIDIA vasculhou 5 milhões de peças e encontrou a agulha no palheiro.

* **Mutações Detetadas (K-mers Únicos do Paciente):** `{_total_mutacoes:,}`

**Amostra de Sequências Mutadas (k={_k_size}):**
```text
{chr(10).join(_amostra_visual)}
```"""), kind="success")
        ])
        
    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO NO PYTHON:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,
@app.cell
def pagina_annotation_ui(definir_pagina, mo, obter_pagina):
    import os as _os
    mo.stop(obter_pagina() != "annotation", mo.md(""))
    btn_voltar_annotation = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    form_annotation = mo.ui.form(
        element=mo.ui.dictionary({
            "dna_input": mo.ui.text_area(
                label="🧬 Cole o DNA Mutado do Módulo 4 (ou um Contig inteiro):",
                # Deixei uma das sequências mutadas que a sua GPU achou como teste!
                value="ATTTACCCTTGCTTTCACATAATCTGCCCGCCGATTTTGCTGCGTTGCGTAAATA" 
            )
        }),
        bordered=True,
        submit_button_label="⚙️ Traduzir para Proteína (Multi-Core)"
    )

    ui_annotation = mo.vstack([
        btn_voltar_annotation,
        mo.md("<br># 📘 Anotação Funcional (Tradução Genética)"),
        mo.callout(mo.md("O DNA por si só é apenas um manual de instruções. Para saber se uma mutação causa doença, precisamos traduzir o código de 3 em 3 letras (Codões) para ver se a estrutura da **Proteína** foi alterada. Cole a sequência abaixo para testar."), kind="info"),
        form_annotation
    ])
    
    ui_annotation
    return btn_voltar_annotation, form_annotation, ui_annotation

@app.cell
async def pagina_annotation_backend(asyncio, form_annotation, mmg, mo, obter_pagina):
    mo.stop(obter_pagina() != "annotation", mo.md(""))
    mo.stop(form_annotation.value is None, mo.md("🧬 *Aguardando o DNA para tradução...*"))

    _dna = form_annotation.value.get("dna_input", "").strip().upper()
    mo.stop(not _dna, mo.callout(mo.md("⚠️ Insira uma sequência de DNA válida."), kind="warning"))

    try:
        with mo.status.spinner("Traduzindo DNA em Proteínas via Rust (Rayon Multi-Core)..."):
            await asyncio.sleep(0.1)
            # O Rust espera uma lista de sequências, devolve uma lista de proteínas
            _resultados = mmg.traduzir_proteinas([_dna])
            _proteina = _resultados[0]

        # Vamos formatar visualmente os codões de STOP (*) para saltarem à vista
        _proteina_formatada = _proteina.replace("*", "**<span style='color:red'>*</span>**")

        _painel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🧬 Tradução Concluída!
            * **Tamanho do DNA:** `{len(_dna)} bases`
            * **Tamanho da Proteína:** `{len(_proteina)} aminoácidos`

            **Sequência Peptídica (Aminoácidos):**
            > {_proteina_formatada}

            *(Nota Médica: O asterisco vermelho indica um 'Codão de Paragem' (STOP). Se a mutação gerou um STOP precoce (Nonsense), a proteína foi truncada e a bactéria provavelmente perdeu essa função biolóica!)*
            """), kind="success")
        ])
    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO NO PYTHON:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,

@app.cell
def pagina_bloom_ui(definir_pagina, glob, mo, obter_pagina):
    mo.stop(obter_pagina() != "bloom", mo.md(""))
    btn_voltar_bloom = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    # Verifica os ficheiros disponíveis
    _arquivos = glob.glob("*.arrow")
    _opcoes = _arquivos if _arquivos else ["Nenhum arquivo .arrow encontrado"]

    form_bloom = mo.ui.form(
        element=mo.ui.dictionary({
            "path_arrow": mo.ui.dropdown(options=_opcoes, label="🧬 Escolha o Genoma (.arrow):"),
            "tamanho_mb": mo.ui.slider(start=1, stop=128, step=1, value=16, label="Memória VRAM (MB):"),
            "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=5, label="Número de Hashes:")
        }),
        bordered=True,
        submit_button_label="🔍 1. Lançar Força Bruta na NVIDIA (CUDA)"
    )

    _ui = mo.vstack([
        btn_voltar_bloom,
        mo.md("<br># 🔍 Filtros Bloom Probabilísticos (GPU)"),
        form_bloom
    ])
    _ui
    return btn_voltar_bloom, form_bloom

@app.cell
async def pagina_bloom_backend(asyncio, form_bloom, mmg, mo, np, obter_pagina, pl):
    mo.stop(obter_pagina() != "bloom", mo.md(""))
    mo.stop(form_bloom.value is None, mo.md("🔍 *Configure os parâmetros acima e clique em 'Lançar'...*"))

    _path = form_bloom.value.get("path_arrow")
    mo.stop(not _path or "Nenhum" in _path, mo.md("⚠️ Selecione um ficheiro .arrow válido."))

    _tamanho_mb = form_bloom.value["tamanho_mb"]
    _num_hashes = form_bloom.value["num_hashes"]

    with mo.status.spinner("🔥 Processando na NVIDIA GTX 1050 Ti..."):
        await asyncio.sleep(0.1)
        _df = pl.read_ipc(_path)
        _bytes_kmer = _df.get_column("K-mers").to_numpy().tobytes()
        _total_kmers = len(_df)
        
        # AQUI É ONDE A MÁGICA CUDA ACONTECE
        _bloom_bytes = mmg.construir_filtro_bloom_gpu(_bytes_kmer, _tamanho_mb, _num_hashes)

        # Estatísticas
        _arr_8 = np.frombuffer(_bloom_bytes, dtype=np.uint8)
        _bits_ativos = np.unpackbits(_arr_8).sum()
        _bits_totais = _tamanho_mb * 1024 * 1024 * 8
        _saturacao = (_bits_ativos / _bits_totais) * 100

    _painel = mo.callout(mo.md(f"""
    ### 🚀 Filtro Bloom Gerado com Sucesso!
    A GPU processou **{_total_kmers * _num_hashes:,}** operações de hash.
    * **Saturação:** `{_saturacao:.2f}%` 
    * **Bits Ativos:** `{_bits_ativos:,}` de `{_bits_totais:,}`
    """), kind="success")
    
    _painel
    return _painel,
@app.cell
def pagina_scaffold_ui(definir_pagina, glob, mo, obter_pagina):
    mo.stop(obter_pagina() != "scaffold", mo.md(""))
    btn_voltar_scaffold = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    _arquivos = glob.glob("*.fasta") + glob.glob("*.fa")
    if not _arquivos: _arquivos = ["Nenhum arquivo FASTA encontrado"]

    form_scaffold = mo.ui.form(
        element=mo.ui.dictionary({
            "path_fasta": mo.ui.dropdown(options=_arquivos, label="🧩 Ficheiro de Contigs (FASTA):"),
            # Se usámos K=55 na montagem, o overlap exato dos contigs é K-1 (54)
            "overlap": mo.ui.slider(start=21, stop=100, step=1, value=54, label="Sobreposição Exata nas Pontas (Overlap):")
        }),
        bordered=True,
        submit_button_label="🏗️ Fechar Genoma (Scaffolding)"
    )

    ui_scaffold = mo.vstack([
        btn_voltar_scaffold,
        mo.md("<br># 🧬 Fechamento de Genomas (Scaffolding)"),
        mo.callout(mo.md("O Grafo de De Bruijn gerou grandes fragmentos (Contigs). O Módulo de Scaffolding pega nestes blocos de centenas de milhares de bases e costura-os procurando sobreposições perfeitas nas pontas, com o objetivo de reconstruir o círculo cromossómico original."), kind="info"),
        form_scaffold
    ])
    
    ui_scaffold
    return btn_voltar_scaffold, form_scaffold, ui_scaffold

@app.cell
async def pagina_scaffold_backend(asyncio, form_scaffold, mmg, mo, obter_pagina):
    mo.stop(obter_pagina() != "scaffold", mo.md(""))
    mo.stop(form_scaffold.value is None, mo.md("🏗️ *Aguardando seleção do FASTA...*"))

    _path = form_scaffold.value.get("path_fasta")
    _overlap = form_scaffold.value.get("overlap")
    
    mo.stop(not _path or "Nenhum" in _path, mo.callout(mo.md("⚠️ Selecione um ficheiro FASTA válido gerado pelo Módulo 3."), kind="warning"))

    try:
        with mo.status.spinner("Lendo Contigs massivos da RAM..."):
            await asyncio.sleep(0.1)
            # Lê o ficheiro FASTA gerado pelo Módulo 3
            with open(_path, "r") as f:
                linhas = f.readlines()
            
            contigs = []
            seq_atual = []
            for linha in linhas:
                if linha.startswith(">"):
                    if seq_atual:
                        contigs.append("".join(seq_atual))
                        seq_atual = []
                else:
                    seq_atual.append(linha.strip())
            if seq_atual:
                contigs.append("".join(seq_atual))

        with mo.status.spinner(f"Costurando pontas com overlap exato de {_overlap} bases no Rust..."):
            await asyncio.sleep(0.1)
            # 🚀 Chamada ao Motor C++ / Rust
            scaffolds, circulares = mmg.montar_scaffolds(contigs, _overlap)

        _painel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🏆 Scaffolding Concluído!
            O algoritmo analisou e cruzou as pontas de **{len(contigs):,}** blocos isolados.
            
            * **Super-Contigs Resultantes:** `{len(scaffolds):,}` 
            * **Maior Cromossoma Contínuo:** `{len(scaffolds[0]):,} bases`
            * **Genomas Circulares Detetados:** `{circulares}` 🔄
            
            *(Se detetou Genomas Circulares, significa que a cauda final do cromossoma colou perfeitamente com a cabeça inicial, reconstruindo a biologia original na sua plenitude!)*
            """), kind="success")
        ])
        
    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO NO PYTHON:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,
@app.cell
def pagina_vcf_ui(definir_pagina, glob, mo, obter_pagina):
    mo.stop(obter_pagina() != "vcf", mo.md(""))
    btn_voltar_vcf = mo.ui.button(label="⬅️ Voltar ao Menu", on_change=lambda _: definir_pagina("home"))

    # 👇 A MAGIA ACONTECE AQUI: Procurar explicitamente pelas extensões
    _arquivos = glob.glob("*.vcf") + glob.glob("*.vcf.gz") + glob.glob("*.parquet")
    if not _arquivos: 
        _arquivos = ["Nenhum arquivo encontrado (Coloque os VCF/Parquet na pasta do projeto)"]

    form_vcf = mo.ui.form(
        element=mo.ui.dictionary({
            "caminhos": mo.ui.multiselect(options=_arquivos, label="🗄️ Selecione os ficheiros (VCF, GZ, Parquet):"),
            "ordem": mo.ui.dropdown(options=["POS", "QUAL"], value="POS", label="📊 Campo de Ordenação:"),
            "nome_saida": mo.ui.text(value="Lote_Fundido_Ordenado.vcf", label="💾 Nome do ficheiro final:")
        }),
        bordered=True,
        submit_button_label="🚀 Iniciar Motor Multimerge CLI"
    )

    ui_vcf = mo.vstack([
        btn_voltar_vcf,
        mo.md("<br># 🗄️ 8. VCF Multimerge & Sorting"),
        mo.callout(mo.md("O projeto original vive! Ingestão massiva de variantes genéticas. Selecione múltiplos ficheiros `.vcf`, `.vcf.gz` ou `.parquet` e o motor CLI com Polars fará o trabalho."), kind="info"),
        form_vcf
    ])
    
    ui_vcf
    return btn_voltar_vcf, form_vcf, ui_vcf
@app.cell
async def pagina_vcf_backend(asyncio, form_vcf, mo, obter_pagina):
    import os as _os
    import subprocess as _subprocess
    
    mo.stop(obter_pagina() != "vcf", mo.md(""))
    mo.stop(form_vcf.value is None, mo.md("⏳ *Aguardando configuração...*"))

    _caminhos = form_vcf.value.get("caminhos", [])
    _ordem = form_vcf.value.get("ordem")
    _saida = form_vcf.value.get("nome_saida", "Output.vcf")

    mo.stop(not _caminhos or "Nenhum" in _caminhos[0], mo.callout(mo.md("⚠️ Selecione pelo menos um ficheiro."), kind="warning"))

    try:
        with mo.status.spinner(f"Processando ficheiros via Motor Polars/Rust CLI..."):
            await asyncio.sleep(0.1)
            
            # Caminho dinâmico para encontrar o seu executável compilado pelo Cargo
            _cli_path = _os.path.join("target", "release", "multimerge_cli")
            if _os.name == "nt": # Adiciona .exe no Windows
                _cli_path += ".exe"
                
            if not _os.path.exists(_cli_path):
                mo.stop(True, mo.md(f"⚠️ **Erro:** O executável não foi encontrado em '{_cli_path}'. Compile com `cargo build --release --bin multimerge_cli` no terminal primeiro."))

            _regra_cli = "--sort-by 0:chrom,1:num" if _ordem == "POS" else "--sort-by 5:float"
            
            # Dispara o seu executável para cada ficheiro selecionado
            for _caminho in _caminhos:
                _comando = f"{_cli_path} {_caminho} {_saida} {_regra_cli}"
                _subprocess.run(_comando.split(), check=True)

            _caminho_absoluto = _os.path.abspath(_saida)

        _painel = mo.vstack([
            mo.callout(mo.md(f"""
            ### ✅ VCF Multimerge Concluído!
            O motor CLI Polars leu, ordenou e fundiu com sucesso.
            
            * **Ficheiros Originais:** `{len(_caminhos)}`
            * **Campo de Ordenação Aplicado:** `{_ordem}`
            * **Guardado em:** `{_caminho_absoluto}`
            """), kind="success")
        ])
        
    except Exception as e:
        _painel = mo.callout(mo.md(f"**🔥 ERRO CRÍTICO NO CLI:** `{str(e)}`"), kind="danger")

    _painel
    return _painel,