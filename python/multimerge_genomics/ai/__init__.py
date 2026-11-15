"""AI for MMGEN: a catalog of the tools, an assistant that explains a result
or turns a problem into a path through the tools, an agent that runs them with
your approval, and a session record.

The catalog, the session record and the privacy helpers need nothing extra.
The assistant and the agent call Claude through the Anthropic SDK
(`pip install "multimerge-genomics[ai]"`, and an API key in ANTHROPIC_API_KEY
or passed in). Experimental, for teaching and research.
"""
from ._catalog import SCREENS, catalog, catalog_text, page_of, tool_definitions
from .privacy import compact, result_for_ai
from .session import Call, Recorder

__all__ = [
    "SCREENS", "catalog", "catalog_text", "page_of", "tool_definitions", "compact", "result_for_ai", "Call", "Recorder",
    "Assistant", "Answer", "AIUnavailable", "ConsentRequired", "Agent", "check_script", "MODELS", "DEFAULT_MODEL",
]


def __getattr__(name):
    # The assistant and the agent load on first use: importing the catalog
    # never needs the Anthropic SDK.
    if name in ("Assistant", "Answer", "AIUnavailable", "ConsentRequired", "check_script", "MODELS", "DEFAULT_MODEL"):
        from . import assistant

        return getattr(assistant, name)
    if name == "Agent":
        from .agent import Agent

        return Agent
    raise AttributeError(f"module 'multimerge_genomics.ai' has no attribute '{name}'")
