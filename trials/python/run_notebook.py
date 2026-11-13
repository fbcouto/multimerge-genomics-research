# Runs the code cells of a notebook, in order, in a single namespace, and shows
# the output of each. Skips the cells that need the browser (files.download).
# Usage: python run_notebook.py <notebook.ipynb>
import json
import sys
import traceback

nb = json.load(open(sys.argv[1], encoding="utf-8"))
ns = {"__name__": "__main__"}
falhas = 0
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    fonte = "".join(c["source"])
    titulo = fonte.splitlines()[0]
    if "files.download" in fonte:
        print(f"\n--- [{i}] {titulo}\n(skipped: needs the browser)")
        continue
    print(f"\n--- [{i}] {titulo}", flush=True)
    try:
        exec(compile(fonte, f"celula_{i}", "exec"), ns)
    except Exception:
        falhas += 1
        traceback.print_exc()
print(f"\n{falhas} cell(s) with errors")
sys.exit(1 if falhas else 0)
