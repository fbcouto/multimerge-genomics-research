"""MMGEN's tools as an MCP server, for Claude (desktop, claude.ai, Claude Code)
and other MCP clients.

    pip install "multimerge-genomics[mcp]"
    mmgen-mcp --root C:/my/data          (or: python -m multimerge_genomics.mcp_server --root ...)

Every tool that takes JSON arguments is offered, with its documentation and a
JSON schema. Paths are relative to `--root` (default: the current folder), and
no path outside it is accepted, for reading or writing. Human-genome tools are
offered too; their results carry the experimental notice. MCP clients ask the
person before each call. Experimental, for teaching and research.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from . import _native
from .ai import _catalog as catalog, privacy
from .ai.agent import _looks_like_path

INSTRUCTIONS = (
    "MMGEN (MultiMerge Genomics): genomics tools for files of any size - validate, sort, index and query VCF/BED/GFF/SAM/BAM, "
    "interval operations, VCF statistics/filter/normalization, FASTQ quality and trimming, k-mer counting, human genetics "
    "(kinship, ancestry, ClinVar findings, pharmacogenomics via PharmCAT, ROH), RNA-seq, CRISPR screens, metagenomics, "
    "resistance genes, outbreak trees, viral mutations, in-silico PCR, 16S amplicons, tumor mutational burden. "
    "Call mmgen_catalog first to see the screens and tools. Paths are relative to the server's working folder. "
    "Human-genome results are experimental, for teaching and research: not a diagnosis nor advice on any treatment."
)


def _inside(root: str, value: str):
    p = os.path.abspath(value if os.path.isabs(value) else os.path.join(root, value))
    try:
        return p if os.path.commonpath([p, root]) == root else None
    except ValueError:
        return None


def resolve_arguments(root: str, schema: dict, arguments: dict) -> dict:
    """Checks the arguments of a call and turns its paths into absolute paths
    inside `root`. Raises ValueError with the reason when it cannot run."""
    props = schema.get("properties", {})
    out = {}
    for k, v in (arguments or {}).items():
        if k not in props:
            raise ValueError(f"no parameter `{k}`; the parameters are: {', '.join(props)}")
        def fix(s):
            if isinstance(s, str) and _looks_like_path(k, s):
                p = _inside(root, s)
                if p is None:
                    raise ValueError(f"`{k}` = {s!r} is outside the working folder {root}")
                return p
            return s
        out[k] = [fix(s) for s in v] if isinstance(v, list) else fix(v)
    missing = [r for r in schema.get("required", []) if r not in out]
    if missing:
        raise ValueError(f"missing required parameter(s): {', '.join(missing)}")
    return out


def call_tool(root: str, name: str, arguments: dict) -> tuple[bool, str]:
    """Runs one tool: (ok, text). The text is the result as JSON (long lists
    cut), or the error."""
    if name == "mmgen_catalog":
        lines = [f"{s[0]}. {s[2]} ({s[3]}): {s[4]}." for s in catalog.SCREENS]
        lines.append("")
        for t in catalog.catalog(include_bytes=False):
            lines.append(f"mmg.{t['signature']} [screen {t['screen']}]{' [human genome]' if t['human'] else ''}: {t['summary']}")
        return True, "\n".join(lines)
    tools = {t["name"]: t for t in catalog.tool_definitions(human=True)}
    if name not in tools:
        return False, f"Unknown tool `{name}`."
    try:
        args = resolve_arguments(root, tools[name]["input_schema"], arguments)
        result = getattr(_native, name)(**args)
    except Exception as err:
        return False, f"{type(err).__name__}: {err}"
    return True, privacy.result_for_ai(result, max_chars=40000)


def build_server(root: str):
    import anyio
    import mcp.types as types
    from mcp.server.lowlevel import Server

    root = os.path.abspath(root)
    defs = catalog.tool_definitions(human=True)
    tool_list = [types.Tool(name="mmgen_catalog", description="The screens of MMGEN and every tool, one line each. Call this first.",
                            input_schema={"type": "object", "properties": {}, "required": []})]
    tool_list += [types.Tool(name=t["name"], description=t["description"], input_schema=t["input_schema"]) for t in defs]

    async def on_list_tools(ctx, params):
        return types.ListToolsResult(tools=tool_list)

    async def on_call_tool(ctx, params):
        ok, text = await anyio.to_thread.run_sync(call_tool, root, params.name, params.arguments or {})
        return types.CallToolResult(content=[types.TextContent(type="text", text=text)], is_error=not ok)

    return Server("mmgen", version=_version(), instructions=INSTRUCTIONS + f" Working folder: {root}",
                  on_list_tools=on_list_tools, on_call_tool=on_call_tool)


def _version() -> str:
    try:
        return str(_native.engine_info()["version"])
    except Exception:
        return ""


async def serve(root: str) -> None:
    from mcp.server.stdio import stdio_server

    server = build_server(root)
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="mmgen-mcp", description="MMGEN's tools as an MCP server (stdio).")
    ap.add_argument("--root", default=os.getcwd(), help="working folder: every path is relative to it and must stay inside it")
    args = ap.parse_args(argv)
    if not os.path.isdir(args.root):
        sys.exit(f"--root {args.root!r} is not a folder")
    try:
        import anyio
        import mcp  # noqa: F401
    except ImportError:
        sys.exit('The MCP server needs the MCP SDK: pip install "multimerge-genomics[mcp]" (or pip install mcp).')
    anyio.run(serve, args.root)


if __name__ == "__main__":
    main()
