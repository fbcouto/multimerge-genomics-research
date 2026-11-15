"""A record of the tools run in a session: what the AI explains ("the last
result of this screen") and what the person can take home as a Python script
that does the same again."""
from __future__ import annotations

import datetime as _dt
import threading
import time
from dataclasses import dataclass, field

from . import _catalog as catalog

_SKIP = {"engine_info", "parquet_options"}


class Bytes:
    """Stands in for a `bytes` argument or result (k-mers, a Bloom filter):
    the record keeps its size, not the data."""

    def __init__(self, n: int):
        self.n = n

    def __repr__(self):
        return f"<{self.n:,} bytes>"


def _shrink(v):
    if isinstance(v, (bytes, bytearray, memoryview)):
        return Bytes(len(v))
    if isinstance(v, (list, tuple)) and v and all(isinstance(x, (bytes, bytearray)) for x in v):
        return type(v)(Bytes(len(x)) for x in v)
    return v


def _has_bytes(v) -> bool:
    return isinstance(v, Bytes) or (isinstance(v, (list, tuple)) and any(isinstance(x, Bytes) for x in v))


@dataclass
class Call:
    name: str
    args: tuple
    kwargs: dict
    result: object = None
    error: str | None = None
    seconds: float = 0.0
    when: str = field(default_factory=lambda: _dt.datetime.now().strftime("%H:%M:%S"))
    context: str | None = None  # the screen it was run from (page key), when known

    @property
    def page(self):
        return self.context or catalog.page_of(self.name)

    @property
    def needs_bytes(self) -> bool:
        """It took the bytes of an earlier call: the script cannot repeat it alone."""
        return any(_has_bytes(a) for a in self.args) or any(_has_bytes(v) for v in self.kwargs.values())

    def python(self) -> str:
        parts = [repr(a) for a in self.args] + [f"{k}={v!r}" for k, v in self.kwargs.items()]
        return f"mmg.{self.name}({', '.join(parts)})"


class Recorder:
    """Stands in for the `multimerge_genomics` module: every function call goes
    to the engine and is recorded (name, arguments, result, time). `context`,
    when set, is stored with each call: the interface puts the screen there."""

    def __init__(self, module, keep: int = 200):
        self._module = module
        self._keep = keep
        self._lock = threading.Lock()
        self.calls: list[Call] = []
        self.context: str | None = None

    def __getattr__(self, name):
        attr = getattr(self._module, name)
        if not callable(attr) or name.startswith("_") or name in _SKIP:
            return attr

        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            call = Call(name, tuple(_shrink(a) for a in args), {k: _shrink(v) for k, v in kwargs.items()},
                        context=self.context)
            try:
                result = attr(*args, **kwargs)
                call.result = _shrink(result)
                return result
            except Exception as err:
                call.error = f"{type(err).__name__}: {err}"
                raise
            finally:
                call.seconds = time.perf_counter() - t0
                with self._lock:
                    self.calls.append(call)
                    del self.calls[: -self._keep]

        wrapper.__name__ = name
        wrapper.__doc__ = attr.__doc__
        return wrapper

    def last(self, page: str | None = None) -> Call | None:
        """The last call that worked, of one screen (page key) or of any. A call
        whose result is only bytes is passed over when there is a better one."""
        with self._lock:
            calls = [c for c in self.calls if c.error is None and (page is None or c.page == page)]
        for c in reversed(calls):
            if not isinstance(c.result, Bytes):
                return c
        return calls[-1] if calls else None

    def clear(self) -> None:
        with self._lock:
            self.calls.clear()

    def script(self) -> str:
        """The session as a Python script that runs the same calls again."""
        lines = [
            "# MMGEN session, " + _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "# Every tool run in the web interface, in order, as Python calls.",
            "# Calls that took the bytes of an earlier call (the k-mer screens) are left as comments.",
            "import multimerge_genomics as mmg",
            "",
        ]
        with self._lock:
            calls = list(self.calls)
        for c in calls:
            screen = catalog.SCREEN_BY_PAGE.get(c.page)
            where = f"screen {screen[0]}, {screen[2]}" if screen else "general"
            note = "  FAILED: " + c.error if c.error else ("  (takes bytes from an earlier call)" if c.needs_bytes else "")
            lines.append(f"# {c.when}  {where}  ({c.seconds:.1f} s)" + note)
            lines.append(("# " if c.error or c.needs_bytes else "") + f"r = {c.python()}")
            lines.append("")
        return "\n".join(lines)
