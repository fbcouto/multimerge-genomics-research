"""What leaves the machine when the AI is used: never a genome or a file, only
a short summary of a tool's result, with long lists cut and file paths reduced
to their file names (a folder name can carry a person's name)."""
from __future__ import annotations

import json
import math
import os
import re

_SEQUENCE = re.compile(r"[A-Za-z*\-]+")
_DNA = re.compile(r"[ACGTUNacgtun\-]+")

# A string that is one path, spaces allowed: it starts like a path (C:\, /,
# ./, ../) and holds no other path after a space.
_WHOLE = re.compile(r"(?:[A-Za-z]:[\\/]|[\\/]|\.{1,2}[\\/])")
_ANOTHER = re.compile(r"\s(?:[A-Za-z]:[\\/]|[\\/]|\.{1,2}[\\/])")
# Paths inside a longer text (a report that names the file it read): the
# folders go, the file name stays. Windows (C:\..., not the s:// of a URL),
# absolute POSIX (/home/..., not the // of a URL) and relative ones that end
# in a file name with an extension (data/maria/x.vcf; not a genotype 0/1,
# Ts/Tv or a ratio 3.5/4.0).
_IN_TEXT = (
    re.compile(r"(?<![A-Za-z])[A-Za-z]:[\\/](?:[^\\/\s'\"`<>|]+[\\/])*"),
    re.compile(r"(?<![\w.:/\\])/(?:[^/\s'\"`<>]+/)+"),
    re.compile(r"(?<![\w.:/\\-])(?:\.{1,2}[\\/])?(?:[\w.-]*[A-Za-z][\w.-]*[\\/])+(?=[\w-]+\.[A-Za-z]\w{0,7}\b)"),
)


def scrub_paths(s: str) -> str:
    """The text with the folders of every path in it removed."""
    for rx in _IN_TEXT:
        s = rx.sub("", s)
    return s


def _shorten_path(s: str) -> str:
    if "/" not in s and "\\" not in s:
        return s
    if len(s) < 400 and "\n" not in s and _WHOLE.match(s) and not _ANOTHER.search(s):
        return re.split(r"[\\/]", s.rstrip("/\\"))[-1] or s
    return scrub_paths(s)


def compact(obj, max_list: int = 30, max_str: int = 1500, _depth: int = 0):
    """A copy of a result small enough to send: lists cut to `max_list`
    items (with a note of how many more), long strings cut, floats rounded,
    paths reduced to their file names."""
    if _depth > 8:
        return "..."
    if isinstance(obj, dict):
        return {str(k): compact(v, max_list, max_str, _depth + 1) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        items = [compact(v, max_list, max_str, _depth + 1) for v in list(obj)[:max_list]]
        if len(obj) > max_list:
            items.append(f"... and {len(obj) - max_list} more")
        return items
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return str(obj)
        return float(f"{obj:.6g}")
    if isinstance(obj, (bytes, bytearray)):
        return f"<{len(obj)} bytes>"
    if isinstance(obj, str):
        if len(obj) > 100 and _SEQUENCE.fullmatch(obj):
            # a DNA or protein sequence: its length says enough, the letters stay here
            kind = "DNA sequence" if _DNA.fullmatch(obj) else "sequence"
            return f"<{kind} of {len(obj):,} letters>"
        s = _shorten_path(obj)
        return s if len(s) <= max_str else s[:max_str] + f"... ({len(s) - max_str} more characters)"
    if obj is None or isinstance(obj, (bool, int)):
        return obj
    return compact(str(obj), max_list, max_str, _depth + 1)


def result_for_ai(result, max_chars: int = 12000) -> str:
    """The result as JSON text for the model, at most `max_chars` long."""
    text = json.dumps(compact(result), ensure_ascii=False, indent=1)
    if len(text) > max_chars:
        text = json.dumps(compact(result, max_list=8, max_str=400), ensure_ascii=False, indent=1)
    if len(text) > max_chars:
        text = text[:max_chars] + "\n... (cut)"
    return text


def call_for_ai(call: str) -> str:
    """A Python call with its file paths reduced to file names."""
    return re.sub(r"'([^']*[\\/][^']*)'", lambda m: repr(_shorten_path(m.group(1)).replace("\\\\", "\\")), call)
