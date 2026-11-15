# Test of the AI layer without a network or an API key: the catalog, the
# privacy helpers, the session record, the assistant's requests and answers
# (with a fake Claude client), the agent's approval loop and path policy, and
# the MCP server end to end (a real stdio session with the MCP SDK).
# Usage: python trials/python/test_ai.py      (needs: anthropic, mcp)
import asyncio
import json
import os
import sys
import tempfile
from types import SimpleNamespace as NS

import multimerge_genomics as mmg
from multimerge_genomics import ai
from multimerge_genomics.ai import agent as agent_mod
from multimerge_genomics.ai import assistant as asst
from multimerge_genomics.ai import knowledge, privacy

failures, done = [], 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {str(detail)[:400]}"))


T = tempfile.mkdtemp(prefix="mmg_ai_")
vcf = os.path.join(T, "a.vcf")
with open(vcf, "w", newline="\n") as fh:
    fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i in range(50):
        fh.write(f"chr1\t{1000 + 37 * i}\t.\tA\tG\t50\tPASS\t.\n")

print("== the catalog")
cat = ai.catalog()
names = {t["name"] for t in cat}
check("every engine function is in the catalog", names == {n for n in dir(mmg._native) if not n.startswith("_") and callable(getattr(mmg._native, n))},
      sorted(names ^ {n for n in dir(mmg._native) if not n.startswith("_")}))
check("73 functions", len(cat) == 73, len(cat))
check("every function but engine_info has a screen", all(t["screen"] for t in cat if t["name"] != "engine_info"),
      [t["name"] for t in cat if not t["screen"]])
check("12 pass raw bytes", sum(not t["json_callable"] for t in cat) == 12)
sg = next(t for t in cat if t["name"] == "sort_genomic")
check("sort_genomic: required and typed parameters", sg["input_schema"]["required"] == ["paths", "output"]
      and sg["input_schema"]["properties"]["index"] == {"type": "boolean"} and sg["screen"] == 9 and sg["outputs"] == ["output"], sg["input_schema"])
check("human-genome tools are marked", next(t for t in cat if t["name"] == "relatedness")["human"]
      and not next(t for t in cat if t["name"] == "gene_regions")["human"] and not sg["human"])
tools_off = {t["name"] for t in ai.tool_definitions(human=False)}
tools_on = {t["name"] for t in ai.tool_definitions(human=True)}
check("tool definitions: human tools only when allowed", "relatedness" not in tools_off and "relatedness" in tools_on and "sort_genomic" in tools_off)
check("tool definitions: no bytes tool", "read_fasta_batch" not in tools_on and "assemble_debruijn" not in tools_on)
check("every tool definition is a valid Claude tool", all(set(t) == {"name", "description", "input_schema"} and t["input_schema"]["type"] == "object"
                                                          for t in ai.tool_definitions(human=True)))
txt = ai.catalog_text()
check("catalog text: screens and signatures", "26. Outbreaks: Trees, Viral Mutations & PCR" in txt and "## mmg.sort_genomic(paths, output" in txt)
check("system prompt: rules, glossary, catalog, notice", all(s in knowledge.system_text() for s in ("teaching assistant of MMGEN", "Short glossary",
                                                                                                     "mmg.vcf_filter(", "take no responsibility")))

print("== privacy")
c = privacy.compact({"path": "C:/Users/maria/genomes/HG002.vcf.gz", "xs": list(range(100)), "f": 1 / 3, "b": b"abc"})
check("paths reduced to the file name", c["path"] == "HG002.vcf.gz", c["path"])
check("long lists cut, floats rounded, bytes described", len(c["xs"]) == 31 and c["xs"][-1] == "... and 70 more" and c["f"] == 0.333333 and c["b"] == "<3 bytes>", c)
check("result_for_ai stays under its limit", len(privacy.result_for_ai({"x": ["y" * 1000] * 500}, max_chars=5000)) <= 5000 + 20)
check("call_for_ai reduces paths", privacy.call_for_ai("mmg.validate_file('D:/data/maria/x.vcf')") == "mmg.validate_file('x.vcf')",
      privacy.call_for_ai("mmg.validate_file('D:/data/maria/x.vcf')"))

print("== the session record")
rec = ai.Recorder(mmg)
r = rec.validate_file(vcf)
check("a recorded call returns the engine's result", r["records"] == 50 and r["valid"])
try:
    rec.validate_file(os.path.join(T, "missing.vcf"))
except Exception:
    pass
rec.engine_info()
check("calls recorded, failures too, engine_info not", [x.name for x in rec.calls] == ["validate_file", "validate_file"] and rec.calls[1].error)
check("last() of a screen skips failures", rec.last("validate") is rec.calls[0] and rec.last("sort") is None)
s = rec.script()
check("the session as a script", "import multimerge_genomics as mmg" in s and f"r = mmg.validate_file({vcf!r})" in s and "# r = mmg.validate_file" in s, s[:400])
rec.context = "sort"  # the interface sets the screen
r2 = rec.validate_file(vcf)
check("the screen a call ran from wins over the function's own", rec.last("sort") is rec.calls[-1] and rec.last("validate") is rec.calls[0])
rec.context = "assembly"
kb = rec.translate_proteins(["ATGGCC"])
fake_bytes = rec.decode_kmers_128(b"\x00" * 32, 31)
check("bytes arguments are recorded by size only", "<32 bytes>" in rec.calls[-1].python() and rec.calls[-1].needs_bytes, rec.calls[-1].python())
check("last() prefers a result that is not bytes", rec.last("assembly").name in ("translate_proteins", "decode_kmers_128"))
s = rec.script()
check("calls that need earlier bytes are comments in the script", "# r = mmg.decode_kmers_128(<32 bytes>, 31)" in s, s[-400:])
check("paths inside a text lose their folders; genotypes, ratios and URLs stay",
      privacy.compact("Read D:\\data\\maria\\HG002.vcf and /home/maria/x.bed and runs/maria/y.vcf.gz: GT 0/1, Ts/Tv 2.1, "
                      "https://github.com/PharmGKB/PharmCAT")
      == "Read HG002.vcf and x.bed and y.vcf.gz: GT 0/1, Ts/Tv 2.1, https://github.com/PharmGKB/PharmCAT",
      privacy.compact("Read D:\\data\\maria\\HG002.vcf and /home/maria/x.bed and runs/maria/y.vcf.gz: GT 0/1, Ts/Tv 2.1, "
                      "https://github.com/PharmGKB/PharmCAT"))
check("DNA sequences never leave as letters", privacy.compact({"contig": "ACGT" * 100})["contig"] == "<DNA sequence of 400 letters>")


# --- a fake Claude client -----------------------------------------------------
def message(text="", tools=(), stop="end_turn", category=None, fallback=False):
    content = [NS(type="thinking", thinking="", signature="sig")]
    if fallback:
        content.append(NS(type="fallback"))
    if text:
        content.append(NS(type="text", text=text))
    for i, (name, inp) in enumerate(tools):
        content.append(NS(type="tool_use", id=f"tu_{i}_{name}", name=name, input=inp))
    return NS(content=content, stop_reason=stop, stop_details=NS(category=category) if category else None, model="claude-opus-5-5",
              usage=NS(input_tokens=1000, output_tokens=100, cache_read_input_tokens=900, cache_creation_input_tokens=0))


class FakeStream:
    def __init__(self, msg):
        self.msg = msg

    async def __aenter__(self):
        return self

    async def __aexit__(self, *a):
        return False

    @property
    def text_stream(self):
        async def gen():
            text = "".join(b.text for b in self.msg.content if b.type == "text")
            for i in range(0, len(text), 7):
                yield text[i:i + 7]
        return gen()

    async def get_final_message(self):
        return self.msg


class FakeClient:
    def __init__(self, script):
        self.script, self.requests = list(script), []
        outer = self

        class Messages:
            def create(self, **kw):
                outer.requests.append(json.loads(json.dumps(kw, default=lambda o: repr(o))))
                return outer.script.pop(0)

            def stream(self, **kw):
                outer.requests.append(json.loads(json.dumps(kw, default=lambda o: repr(o))))
                return FakeStream(outer.script.pop(0))

        self.beta = NS(messages=Messages())


print("== the assistant")
fake = FakeClient([message("The file has 50 SNPs; Ts/Tv is high because all are A>G.")])
asst.client_factory = lambda key, is_async: fake
a = ai.Assistant()
ans = a.explain(f"mmg.vcf_stats({vcf!r})", mmg.vcf_stats(vcf))
req = fake.requests[0]
check("model, effort and fallbacks by default", req["model"] == "claude-opus-5-5" and req["output_config"] == {"effort": "medium"}
      and req["fallbacks"] == "default" and req["betas"] == ["server-side-fallback-2026-07-01"], {k: req.get(k) for k in ("model", "betas", "fallbacks")})
check("the system prompt is cached", req["system"][0]["cache_control"] == {"type": "ephemeral"} and "MMGEN" in req["system"][0]["text"])
check("no thinking setting sent (adaptive on this model)", "thinking" not in req)
user = req["messages"][0]["content"]
check("only the file name and the result summary go out", "a.vcf" in user and T not in user and '"snps": 50' in user, user[:300])
check("the answer is read", ans.text.startswith("The file has 50 SNPs") and not ans.refused and ans.usage["cache_read_input_tokens"] == 900)
fake2 = FakeClient([message("", stop="refusal", category="bio")])
asst.client_factory = lambda key, is_async: fake2
ans2 = ai.Assistant().explain("mmg.relatedness(['a.vcf'])", {"pairs": []})
check("a refusal is reported with its category", ans2.refused and ans2.refusal_category == "bio")
check("a human-genome call asks for the human rules", "human-genome result" in fake2.requests[0]["messages"][0]["content"])
fake3 = FakeClient([message("ok")])
asst.client_factory = lambda key, is_async: fake3
ai.Assistant(model="claude-haiku-5-5").ask("hi")
check("Haiku: no server-side fallbacks", "fallbacks" not in fake3.requests[0] and "betas" not in fake3.requests[0])
plan_text = ("1. **The problem** ...\n```python\nimport multimerge_genomics as mmg\nr = mmg.sort_genomic(['a.vcf'], 'a.sorted.vcf.gz', index=True)\n"
             "s = mmg.vcf_stats('a.sorted.vcf.gz')\n```\n")
fake4 = FakeClient([message(plan_text)])
asst.client_factory = lambda key, is_async: fake4
seen = []


async def grow(t):
    seen.append(t)


p = asyncio.run(ai.Assistant().plan_stream("I want to sort and summarize my VCF", ["a.vcf"], on_text=grow))
check("streaming: the text grows, the final answer is whole", len(seen) > 3 and seen[-1] == plan_text and p.text == plan_text, len(seen))
check("a correct script passes the check", p.warnings == [], p.warnings)
check("the plan prompt lists the files and the sections", "- a.vcf" in fake4.requests[0]["messages"][0]["content"] and "**The script**" in fake4.requests[0]["messages"][0]["content"])
bad = "```python\nmmg.sort_genomics(['a'])\nmmg.vcf_stats('a', colour=1)\nx = (\n```"
w = ai.check_script(bad.replace("x = (\n", "")) + ai.check_script("```python\nx = (\n```")
check("an invented function, parameter or broken script is flagged", any("sort_genomics" in x for x in w) and any("colour" in x for x in w)
      and any("does not parse" in x for x in w), w)

print("== the agent")
root = os.path.join(T, "work")
os.makedirs(root)
os.replace(vcf, os.path.join(root, "a.vcf"))
fake5 = FakeClient([
    message("I will validate the file and read one outside.", tools=[("validate_file", {"path": "a.vcf"}),
                                                                      ("vcf_stats", {"path": "../secret.vcf"}),
                                                                      ("relatedness", {"vcfs": ["a.vcf"]})], stop="tool_use"),
    message("Validated: 50 records, valid. Done.", stop="end_turn"),
])
asst.client_factory = lambda key, is_async: fake5
ag = ai.Agent("Check my VCF", root)
check("the first message lists the files", "- a.vcf" in ag.messages[0]["content"])
r1 = ag.step()
probs = [p.problem for p in ag.pending]
check("three proposals: one fine, one outside the folder, one human tool (off)", len(ag.pending) == 3 and probs[0] is None
      and "outside the working folder" in probs[1] and "not available" in probs[2], probs)
check("the agent sends the tools without the human ones", "relatedness" not in {t["name"] for t in fake5.requests[0]["tools"]})
res = ag.run_pending()
check("only the allowed call runs", [ok for _, ok, _ in res] == [True, False, False] and '"records": 50' in res[0][2], res)
last_user = ag.messages[-1]["content"]
check("every tool_use gets its tool_result, in one message", [b["tool_use_id"] for b in last_user] == ["tu_0_validate_file", "tu_1_vcf_stats", "tu_2_relatedness"]
      and [b["is_error"] for b in last_user] == [False, True, True])
check("the assistant's content goes back unchanged (thinking included)", ag.messages[1]["content"][0].type == "thinking")
ag.step()
check("the final answer ends the run", ag.done and ag.final_text.startswith("Validated") and not ag.pending)
check("the transcript", "Round 1" in ag.transcript() and "Final answer" in ag.transcript())
fake6 = FakeClient([message("", tools=[("sort_genomic", {"paths": ["a.vcf"], "output": "a.sorted.vcf"})], stop="tool_use"), message("Fine, I stop.")])
asst.client_factory = lambda key, is_async: fake6
ag2 = ai.Agent("Sort it", root)
ag2.step()
ag2.decline()
check("a declined round is reported to the model as declined", ag2.messages[-1]["content"][0]["is_error"] and "declined" in ag2.messages[-1]["content"][0]["content"]
      and not os.path.exists(os.path.join(root, "a.sorted.vcf")))
ag3 = ai.Agent("x", root, allow_human=True)
check("human tools when allowed", "relatedness" in ag3.tool_names)
check("region strings and sample names are not paths", not agent_mod._looks_like_path("regions", "chr1:100-200")
      and not agent_mod._looks_like_path("samples", "HG002") and agent_mod._looks_like_path("vcf", "x.vcf"))

print("== the MCP server")
from multimerge_genomics import mcp_server as ms
ok, text = ms.call_tool(root, "vcf_stats", {"path": "a.vcf"})
check("call_tool runs a tool inside the folder", ok and '"records": 50' in text, text[:200])
ok2, text2 = ms.call_tool(root, "vcf_stats", {"path": "../a.vcf"})
check("and refuses a path outside it", not ok2 and "outside the working folder" in text2, text2)
ok3, text3 = ms.call_tool(root, "vcf_stats", {"pth": "a.vcf"})
check("and an unknown parameter", not ok3 and "no parameter `pth`" in text3, text3)
ok4, text4 = ms.call_tool(root, "mmgen_catalog", {})
check("the catalog tool", ok4 and "1. QC & Ingestion" in text4 and "mmg.sort_genomic(" in text4)


async def mcp_session():
    from mcp.client.session import ClientSession
    from mcp.client.stdio import StdioServerParameters, stdio_client

    params = StdioServerParameters(command=sys.executable, args=["-m", "multimerge_genomics.mcp_server", "--root", root],
                                   env={k: v for k, v in os.environ.items() if k != "SSLKEYLOGFILE"})
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            init = await session.initialize()
            tools = await session.list_tools()
            good = await session.call_tool("validate_file", {"path": "a.vcf"})
            outside = await session.call_tool("validate_file", {"path": "../a.vcf"})
            return init, tools, good, outside


try:
    init, tools, good, outside = asyncio.run(asyncio.wait_for(mcp_session(), 120))
    n = len(ai.tool_definitions(human=True)) + 1
    check("MCP: the server introduces itself", init.server_info.name == "mmgen" and "Working folder" in (init.instructions or ""))
    check(f"MCP: {n} tools listed", len(tools.tools) == n and any(t.name == "mmgen_catalog" for t in tools.tools), len(tools.tools))
    check("MCP: a tool call works over stdio", not good.is_error and '"records": 50' in good.content[0].text, good.content[0].text[:200])
    check("MCP: a path outside the folder is an error", outside.is_error and "outside" in outside.content[0].text)
except Exception as err:
    check("MCP session over stdio", False, repr(err))

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
