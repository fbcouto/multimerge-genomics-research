"""The AI assistant: explain a result, or turn a problem described in plain
words into a path through MMGEN's tools (with a Python script).

It calls Claude through the official Anthropic SDK (`pip install
"multimerge-genomics[ai]"`). The key comes from `api_key`, or, when that is
empty, from the environment (`ANTHROPIC_API_KEY`). Only what `privacy` lets
out is sent: the screen, the call with file names, a short summary of the
result and the person's question; never a genome or a file.
"""
from __future__ import annotations

import ast
import time
from dataclasses import dataclass, field

from . import _catalog as catalog, knowledge, privacy

DEFAULT_MODEL = "claude-opus-5-5"
MODELS = {
    "Claude Opus 5.5 (default)": "claude-opus-5-5",
    "Claude Sonnet 5.5": "claude-sonnet-5-5",
    "Claude Haiku 5.5": "claude-haiku-5-5",
}
# If the model's safety classifiers decline a request, the API can run it on
# another model inside the same call ("fallbacks"); Haiku has no such fallback.
FALLBACK_BETA = "server-side-fallback-2026-07-01"
_WITH_FALLBACKS = {"claude-opus-5-5", "claude-sonnet-5-5", "claude-fable-5-2", "claude-fable-5-1", "claude-opus-5"}

# Tests (and offline demonstrations) replace this: client_factory(api_key, is_async) -> client.
client_factory = None


class AIUnavailable(RuntimeError):
    """The AI cannot be used here: the SDK is missing, or there is no key."""


class ConsentRequired(RuntimeError):
    """The person has not agreed to send result summaries to the AI provider."""


@dataclass
class Answer:
    text: str
    refused: bool = False
    refusal_category: str | None = None
    stop_reason: str | None = None
    model: str | None = None
    fallback_used: bool = False
    usage: dict = field(default_factory=dict)
    seconds: float = 0.0
    warnings: list = field(default_factory=list)


def _client(api_key: str | None, is_async: bool):
    if client_factory is not None:
        return client_factory(api_key, is_async)
    try:
        import anthropic
    except ImportError as err:
        raise AIUnavailable('The AI needs the Anthropic SDK: pip install "multimerge-genomics[ai]" (or pip install anthropic).') from err
    try:
        return anthropic.AsyncAnthropic(api_key=api_key or None) if is_async else anthropic.Anthropic(api_key=api_key or None)
    except Exception as err:  # no key anywhere
        raise AIUnavailable(
            "No API key: set ANTHROPIC_API_KEY, or type a key in the AI settings (screen 29). Keys come from "
            "console.anthropic.com; using the API is paid."
        ) from err


def _answer_from(message, seconds: float) -> Answer:
    text = "".join(getattr(b, "text", "") for b in message.content if getattr(b, "type", "") == "text")
    stop = getattr(message, "stop_reason", None)
    details = getattr(message, "stop_details", None)
    usage = getattr(message, "usage", None)
    u = {}
    for k in ("input_tokens", "output_tokens", "cache_read_input_tokens", "cache_creation_input_tokens"):
        v = getattr(usage, k, None) if usage is not None else None
        if v is not None:
            u[k] = v
    return Answer(
        text=text,
        refused=stop == "refusal",
        refusal_category=getattr(details, "category", None) if details is not None else None,
        stop_reason=stop,
        model=getattr(message, "model", None),
        fallback_used=any(getattr(b, "type", "") == "fallback" for b in message.content),
        usage=u,
        seconds=seconds,
    )


class Assistant:
    def __init__(self, api_key: str | None = None, model: str = DEFAULT_MODEL, effort: str = "medium",
                 max_tokens: int = 16000):
        self.api_key = api_key or None
        self.model = model or DEFAULT_MODEL
        self.effort = effort
        self.max_tokens = max_tokens

    # -- the request ---------------------------------------------------------
    def _request(self, user_text: str, max_tokens: int | None = None, **extra) -> dict:
        kw = dict(
            model=self.model,
            max_tokens=max_tokens or self.max_tokens,
            # the rules, the glossary, the catalog and the README: the same on
            # every request, so cached
            system=[{"type": "text", "text": knowledge.system_text(), "cache_control": {"type": "ephemeral"}}],
            messages=[{"role": "user", "content": user_text}],
            output_config={"effort": self.effort},
        )
        if self.model in _WITH_FALLBACKS:
            kw["betas"] = [FALLBACK_BETA]
            kw["fallbacks"] = "default"
        kw.update(extra)
        return kw

    def ask(self, user_text: str) -> Answer:
        client = _client(self.api_key, False)
        t0 = time.perf_counter()
        message = client.beta.messages.create(**self._request(user_text))
        return _answer_from(message, time.perf_counter() - t0)

    async def ask_stream(self, user_text: str, on_text=None) -> Answer:
        """Streams the answer; `on_text(text_so_far)` is called as it grows."""
        client = _client(self.api_key, True)
        t0 = time.perf_counter()
        acc = ""
        async with client.beta.messages.stream(**self._request(user_text)) as stream:
            async for piece in stream.text_stream:
                acc += piece
                if on_text is not None:
                    await on_text(acc)
            message = await stream.get_final_message()
        return _answer_from(message, time.perf_counter() - t0)

    # -- 1. explain a result -------------------------------------------------
    @staticmethod
    def explain_prompt(call: str, result, question: str = "", level: str = "student") -> str:
        name = call.split("(")[0].replace("mmg.", "").strip()
        page = catalog.page_of(name)
        screen = catalog.SCREEN_BY_PAGE.get(page)
        where = f"screen {screen[0]}, {screen[2]}" if screen else "the Python API"
        human = any(t["name"] == name and t["human"] for t in catalog.catalog())
        parts = [
            f"I ran a tool of MMGEN ({where}). Level: {level}.",
            "The call (file paths reduced to file names):",
            "```python\n" + privacy.call_for_ai(call) + "\n```",
            "The result (a summary; long lists are cut):",
            "```json\n" + privacy.result_for_ai(result) + "\n```",
        ]
        if question.strip():
            parts.append("My question: " + question.strip())
        else:
            parts.append(
                "Explain this result: what the tool did, what the main numbers mean, whether they look as expected, "
                "what could be wrong, and what I could do next with MMGEN."
            )
        if human:
            parts.append("This is a human-genome result: keep to the rules for human genomes.")
        return "\n\n".join(parts)

    def explain(self, call: str, result, question: str = "", level: str = "student") -> Answer:
        return self.ask(self.explain_prompt(call, result, question, level))

    async def explain_stream(self, call: str, result, question: str = "", level: str = "student", on_text=None) -> Answer:
        return await self.ask_stream(self.explain_prompt(call, result, question, level), on_text)

    # -- 2. describe your problem ---------------------------------------------
    @staticmethod
    def plan_prompt(problem: str, files=(), level: str = "student") -> str:
        listed = "\n".join(f"- {f}" for f in list(files)[:200]) or "(no file listed)"
        return "\n\n".join([
            f"Level: {level}.",
            "My problem, in my own words:\n" + problem.strip(),
            "The files I have (names only):\n" + listed,
            "Please answer with these sections:\n"
            "1. **The problem**: what I am really asking, explained in detail (the biology, the data, what a good answer looks like).\n"
            "2. **The path**: numbered steps through MMGEN, each with the screen and the function, the inputs, the parameters "
            "that matter and what to expect.\n"
            "3. **Watch out**: what could go wrong, and what MMGEN cannot do here (and what kind of tool could).\n"
            "4. **The script**: one Python block that does the path with `import multimerge_genomics as mmg`, using only "
            "functions and parameters from the catalog, with my file names.",
        ])

    def plan(self, problem: str, files=(), level: str = "student") -> Answer:
        ans = self.ask(self.plan_prompt(problem, files, level))
        ans.warnings = check_script(ans.text)
        return ans

    async def plan_stream(self, problem: str, files=(), level: str = "student", on_text=None) -> Answer:
        ans = await self.ask_stream(self.plan_prompt(problem, files, level), on_text)
        ans.warnings = check_script(ans.text)
        return ans


def python_blocks(text: str) -> list[str]:
    blocks, inside, cur = [], False, []
    for line in text.splitlines():
        if line.strip().startswith("```"):
            if inside:
                blocks.append("\n".join(cur))
                cur, inside = [], False
            elif line.strip().lower().startswith(("```python", "```py")):
                inside = True
            continue
        if inside:
            cur.append(line)
    return blocks


def check_script(text: str) -> list[str]:
    """Checks the Python blocks of an answer against the catalog: every
    `mmg.<function>` must exist and take the keyword arguments it is given.
    Returns the problems found (empty when all is well)."""
    known = {t["name"]: t for t in catalog.catalog()}
    problems = []
    for code in python_blocks(text):
        try:
            tree = ast.parse(code)
        except SyntaxError as err:
            problems.append(f"The script does not parse (line {err.lineno}): {err.msg}.")
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) \
                    and node.func.value.id == "mmg":
                name = node.func.attr
                if name not in known:
                    problems.append(f"`mmg.{name}` is not a function of MMGEN.")
                    continue
                params = set(known[name]["input_schema"]["properties"])
                for kw in node.keywords:
                    if kw.arg is not None and kw.arg not in params:
                        problems.append(f"`mmg.{name}` has no parameter `{kw.arg}`.")
    return problems
