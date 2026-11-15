"""The agent: the AI runs MMGEN's tools on a folder of files, round by round,
and nothing runs until the person approves it.

Each round, the model proposes tool calls; `pending` holds them. The person
approves (`run_pending`) or declines (`decline`); the results go back to the
model, which proposes the next round or writes its final answer. Every path
an argument names must lie inside the working folder (`root`), which is also
where outputs are written; human-genome tools are left out unless allowed.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field

from .. import _native
from . import _catalog as catalog, privacy
from .assistant import Assistant, _answer_from, _client

AGENT_RULES = """You are working as an agent on a folder of files, with MMGEN's tools.
- Plan briefly, then call the tools. The person approves each round of calls before it runs, so propose one coherent
  round at a time (several calls in one round are fine when they do not depend on each other).
- Paths are relative to the working folder; write outputs there too, with clear names. Never overwrite an input.
- After the results come back, check them (records read, lines skipped, warnings) before going on.
- When done, write a short report: what was run, the main numbers, what they mean, what to check, and the Python calls
  that reproduce it.
"""


@dataclass
class Proposal:
    id: str
    name: str
    input: dict
    problem: str | None = None  # why it cannot run (policy), shown before approval


@dataclass
class Round:
    text: str
    proposals: list = field(default_factory=list)
    results: list = field(default_factory=list)  # (name, ok, summary)


class Agent:
    def __init__(self, task: str, root: str, assistant: Assistant | None = None, allow_human: bool = False,
                 max_rounds: int = 8, engine=None):
        self.task = task
        # what runs the tools: the engine itself, or the interface's session
        # record (which then lists the agent's calls too)
        self.engine = engine if engine is not None else _native
        self.root = os.path.abspath(root)
        self.assistant = assistant or Assistant()
        self.allow_human = allow_human
        self.max_rounds = max_rounds
        self.tools = catalog.tool_definitions(human=allow_human)
        self.tool_names = {t["name"] for t in self.tools}
        self.messages: list = [{"role": "user", "content": self._first_message()}]
        self.rounds: list[Round] = []
        self.pending: list[Proposal] = []
        self.done = False
        self.final_text = ""
        self.refused = False
        self.usage = {"input_tokens": 0, "output_tokens": 0}

    def _first_message(self) -> str:
        files = []
        for base, dirs, names in os.walk(self.root):
            dirs[:] = [d for d in dirs if not d.startswith((".", "_")) and d not in ("target", "__pycache__")]
            for n in names:
                rel = os.path.relpath(os.path.join(base, n), self.root).replace(os.sep, "/")
                files.append(rel)
                if len(files) >= 300:
                    break
            if len(files) >= 300 or base.count(os.sep) - self.root.count(os.sep) >= 3:
                dirs[:] = []
        listing = "\n".join(f"- {f}" for f in sorted(files)) or "(empty)"
        return f"{AGENT_RULES}\nThe task:\n{self.task.strip()}\n\nFiles in the working folder:\n{listing}"

    # -- one call to the model ------------------------------------------------
    def step(self) -> Round:
        """Asks the model for the next round. Fills `pending`, or ends the run."""
        if len(self.rounds) >= self.max_rounds:
            self.done = True
            self.final_text = f"Stopped after {self.max_rounds} rounds (the limit)."
            return Round(self.final_text)
        client = _client(self.assistant.api_key, False)
        t0 = time.perf_counter()
        req = self.assistant._request("", tools=self.tools)
        req["messages"] = self.messages
        message = client.beta.messages.create(**req)
        ans = _answer_from(message, time.perf_counter() - t0)
        for k in self.usage:
            self.usage[k] += ans.usage.get(k, 0)
        # the whole content goes back unchanged, thinking blocks included
        self.messages.append({"role": "assistant", "content": message.content})
        rnd = Round(ans.text)
        if ans.refused:
            self.refused, self.done = True, True
            self.final_text = ans.text or f"The AI declined this task ({ans.refusal_category or 'no category given'})."
        elif ans.stop_reason == "tool_use":
            for b in message.content:
                if getattr(b, "type", "") == "tool_use":
                    p = Proposal(b.id, b.name, dict(b.input or {}))
                    p.problem = self.check(p)
                    rnd.proposals.append(p)
            self.pending = rnd.proposals
        else:
            self.done = True
            self.final_text = ans.text
        self.rounds.append(rnd)
        return rnd

    # -- the policy -----------------------------------------------------------
    def _inside(self, value: str) -> str | None:
        p = value if os.path.isabs(value) else os.path.join(self.root, value)
        p = os.path.abspath(p)
        try:
            ok = os.path.commonpath([p, self.root]) == self.root
        except ValueError:  # another drive
            ok = False
        return p if ok else None

    def check(self, p: Proposal) -> str | None:
        """Why a proposed call cannot run, or None."""
        if p.name not in self.tool_names:
            return f"`{p.name}` is not available to the agent" + ("" if self.allow_human else " (human-genome tools are off)") + "."
        props = next(t for t in self.tools if t["name"] == p.name)["input_schema"]["properties"]
        for k, v in p.input.items():
            if k not in props:
                return f"`{p.name}` has no parameter `{k}`."
            for s in (v if isinstance(v, list) else [v]):
                if isinstance(s, str) and _looks_like_path(k, s) and self._inside(s) is None:
                    return f"`{k}` = {s!r} is outside the working folder."
        return None

    def _resolve(self, p: Proposal) -> dict:
        out = {}
        for k, v in p.input.items():
            if isinstance(v, list):
                out[k] = [self._inside(s) if isinstance(s, str) and _looks_like_path(k, s) else s for s in v]
            elif isinstance(v, str) and _looks_like_path(k, v):
                out[k] = self._inside(v)
            else:
                out[k] = v
        return out

    # -- the person's decision -------------------------------------------------
    def run_pending(self) -> list:
        """Runs the approved round and sends the results back to the model."""
        results, blocks = [], []
        for p in self.pending:
            if p.problem:
                ok, text = False, "Not run: " + p.problem
            else:
                try:
                    r = getattr(self.engine, p.name)(**self._resolve(p))
                    ok, text = True, privacy.result_for_ai(r)
                except Exception as err:
                    ok, text = False, f"{type(err).__name__}: {err}"
            results.append((p.name, ok, text))
            blocks.append({"type": "tool_result", "tool_use_id": p.id, "content": text, "is_error": not ok})
        if self.rounds:
            self.rounds[-1].results = results
        self.messages.append({"role": "user", "content": blocks})
        self.pending = []
        return results

    def decline(self, reason: str = "The person declined this round.") -> None:
        blocks = [{"type": "tool_result", "tool_use_id": p.id, "content": reason, "is_error": True} for p in self.pending]
        if self.rounds:
            self.rounds[-1].results = [(p.name, False, reason) for p in self.pending]
        self.messages.append({"role": "user", "content": blocks})
        self.pending = []

    def transcript(self) -> str:
        lines = [f"# Task\n{self.task}\n"]
        for i, r in enumerate(self.rounds, 1):
            lines.append(f"## Round {i}\n{r.text}".rstrip())
            for p in r.proposals:
                lines.append(f"- proposed: mmg.{p.name}({json.dumps(p.input, ensure_ascii=False)})" + (f"  [{p.problem}]" if p.problem else ""))
            for name, ok, _ in r.results:
                lines.append(f"- {'ran' if ok else 'not run/failed'}: {name}")
        if self.done:
            lines.append("\n## Final answer\n" + self.final_text)
        return "\n".join(lines)


_PATH_PARAMS = {
    "path", "paths", "input", "input2", "output", "output2", "output_path", "output_dir", "vcf", "vcfs", "bed", "bams", "fastqs",
    "genomes", "inputs", "a", "b", "annotation", "reference", "panel", "model", "table", "taxonomy", "catalog", "catalog_table",
    "count_table", "counts", "library", "genome", "jar", "positions", "clinvar", "scoring_file", "child", "father", "mother",
    "asv_fasta", "feature_table", "tracking", "metadata", "samples_table", "sequences", "tree_output", "matrix_output",
    "amplicons_output", "reads_output", "summary_output", "pcoa_output", "curves_output", "callable_bed", "callable_beds",
    "databases", "samples", "regions", "host_free", "host_free2", "assignments", "details", "species_reference", "outside_calls",
    "allele_frequencies", "temp_dir", "output_prefix", "file",
}


def _looks_like_path(param: str, value: str) -> bool:
    if param in ("regions",) and ":" in value and not os.path.exists(value):
        return False  # "chr1:1000-2000"
    if param == "samples" and not any(value.lower().endswith(e) for e in (".tsv", ".txt", ".csv", ".fastq", ".fq", ".gz", ".files")):
        return False  # sample names, not files
    return param in _PATH_PARAMS or "/" in value or "\\" in value
