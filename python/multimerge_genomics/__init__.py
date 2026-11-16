"""MultiMerge Genomics: tools for genomic files of any size.

Every function of the engine (written in Rust) is here: validate, normalize
chromosome names, sort, index and query by region, interval operations, VCF
statistics, filtering and normalization, FASTQ quality and trimming, k-mer
counting, the human-genetics, research and microbiology tools, and the k-mer
modules with optional CUDA acceleration. Start with `engine_info()`; every
function has a docstring.

Optional parts, written in Python:

- `multimerge_genomics.ai`: a catalog of the tools for AI assistants
  (`ai.catalog()`), and an assistant that explains a result, turns a described
  problem into a path through the tools, or runs the tools as an agent with
  your approval (`pip install "multimerge-genomics[ai]"`).
- `multimerge_genomics.mcp_server`: the tools as an MCP server, for Claude and
  other MCP clients (`pip install "multimerge-genomics[mcp]"`, then `mmgen-mcp`).
"""
from . import _native
from ._native import *  # noqa: F401,F403  (the engine's functions)

__all__ = [name for name in dir(_native) if not name.startswith("_")] + ["PGX_SCHEMA_PATH"]

# The JSON Schema every pharmacogenomics result follows (pgx_run, either engine).
import os as _os

PGX_SCHEMA_PATH = _os.path.join(_os.path.dirname(__file__), "pgx_result.schema.json")


def __getattr__(name):
    # `multimerge_genomics.ai` loads on first use, so that importing the tools
    # never needs the AI dependencies.
    if name == "ai":
        import importlib

        return importlib.import_module(".ai", __name__)
    raise AttributeError(f"module 'multimerge_genomics' has no attribute '{name}'")
