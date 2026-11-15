#!/usr/bin/env python3
"""Creates a new notebook from the MMGEN notebook template.

Usage:   python notebooks/new_problem.py <folder>/<name> ["Title of the problem"]
Example: python notebooks/new_problem.py my_problems/coverage_of_my_panel "Coverage of my gene panel"

The folder is created inside notebooks/ when it does not exist; an existing notebook is never overwritten.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "_template" / "new_problem.ipynb"
PART = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def main(argv):
    if len(argv) not in (2, 3) or argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        return 0 if len(argv) == 2 else 2
    rel = argv[1].replace("\\", "/").strip("/")
    if rel.endswith(".ipynb"):
        rel = rel[: -len(".ipynb")]
    parts = rel.split("/")
    if len(parts) < 2 or not all(PART.fullmatch(p) for p in parts):
        print(f"Give <folder>/<name> with letters, digits, '_', '-' or '.', e.g. my_problems/coverage_of_my_panel (got {argv[1]!r}).")
        return 2
    target = HERE.joinpath(*parts).with_suffix(".ipynb")
    if target.exists():
        print(f"{target} already exists: choose another name (nothing was changed).")
        return 1
    title = argv[2].strip() if len(argv) == 3 and argv[2].strip() else parts[-1].replace("_", " ").replace("-", " ").capitalize()

    nb = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    slug = re.sub(r"[^A-Za-z0-9-]", "-", parts[-1])[:50]
    for i, cell in enumerate(nb["cells"]):
        text = ("".join(cell["source"])
                .replace("{{TITLE}}", title)
                .replace("{{NOTEBOOK_PATH}}", "/".join(parts) + ".ipynb")
                .replace("{{NAME}}", parts[-1]))
        cell["source"] = text.splitlines(keepends=True)
        cell["id"] = f"{slug}-{i:02d}"
        if cell["cell_type"] == "code":
            cell["outputs"], cell["execution_count"] = [], None
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    print(f"Created {target}")
    print("Open it in Jupyter or Colab, replace the texts in italics, and keep the sections.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
