# Evaluation of the AI assistant against real Claude: the cases of
# trials/ai/cases.json (problems in plain words, and results to explain), each
# answer checked automatically (the functions it should name, a valid script,
# the notice on human genomes, no invented function).
#
# It calls the Claude API, which is PAID: it needs ANTHROPIC_API_KEY and runs
# only with --yes. Without --yes it lists the cases and stops.
#
# Usage: python trials/python/ai_eval.py [--yes] [--model claude-opus-5-5] [--effort medium]
#                                        [--only names,paternity] [--out answers.json]
import argparse
import json
import os
import sys
import time

from multimerge_genomics import ai

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(HERE, "..", "ai", "cases.json")


def grade(case, text):
    """The problems found in one answer (empty when it passes)."""
    low = text.lower()
    problems = []
    for group in case.get("functions", []):
        if not any(f in text for f in group):
            problems.append("names none of " + " / ".join(group))
    for ph in case.get("phrases", []):
        if ph.lower() not in low:
            problems.append(f"lacks {ph!r}")
    anyp = case.get("any_phrases", [])
    if anyp and not any(ph.lower() in low for ph in anyp):
        problems.append("lacks all of " + " / ".join(repr(p) for p in anyp))
    for ph in case.get("forbidden", []):
        if ph.lower() in low:
            problems.append(f"says {ph!r}")
    problems += ai.check_script(text)
    return problems


def main():
    ap = argparse.ArgumentParser(description="Evaluate the AI assistant against real Claude (paid API calls).")
    ap.add_argument("--yes", action="store_true", help="really call the API (it is paid)")
    ap.add_argument("--model", default=ai.DEFAULT_MODEL)
    ap.add_argument("--effort", default="medium")
    ap.add_argument("--only", default="", help="comma-separated case ids")
    ap.add_argument("--out", default="", help="save every answer to this JSON file")
    args = ap.parse_args()

    data = json.load(open(CASES, encoding="utf-8"))
    only = {x.strip() for x in args.only.split(",") if x.strip()}
    cases = [("plan", c) for c in data["plan"]] + [("explain", c) for c in data["explain"]]
    cases = [(k, c) for k, c in cases if not only or c["id"] in only]
    print(f"{len(cases)} case(s), model {args.model}, effort {args.effort}")
    if not args.yes:
        for kind, c in cases:
            print(f"  {kind:8s} {c['id']}")
        print("Nothing sent. Add --yes to call the API (paid; each case is one request).")
        return 0
    if not os.environ.get("ANTHROPIC_API_KEY"):
        sys.exit("Set ANTHROPIC_API_KEY first.")

    assistant = ai.Assistant(model=args.model, effort=args.effort)
    passed, rows, tokens = 0, [], {"input_tokens": 0, "output_tokens": 0, "cache_read_input_tokens": 0}
    t0 = time.perf_counter()
    for kind, c in cases:
        try:
            if kind == "plan":
                ans = assistant.plan(c["problem"], c.get("files", []), c.get("level", "student"))
            else:
                ans = assistant.explain(c["call"], c["result"], c.get("question", ""), c.get("level", "student"))
        except Exception as err:
            print(f"  ERROR  {kind} {c['id']}: {type(err).__name__}: {err}")
            rows.append({"id": c["id"], "kind": kind, "error": str(err)})
            continue
        for k in tokens:
            tokens[k] += ans.usage.get(k, 0)
        problems = ["refused (" + str(ans.refusal_category) + ")"] if ans.refused else grade(c, ans.text)
        ok = not problems
        passed += ok
        print(f"  {'ok   ' if ok else 'FAIL '}  {kind:8s} {c['id']:14s} {ans.seconds:5.1f} s  "
              f"{ans.usage.get('output_tokens', 0):5d} tokens out" + ("" if ok else "  " + "; ".join(problems)))
        rows.append({"id": c["id"], "kind": kind, "ok": ok, "problems": problems, "text": ans.text, "model": ans.model,
                     "fallback_used": ans.fallback_used, "usage": ans.usage, "seconds": ans.seconds})
    print(f"\n{passed} of {len(cases)} cases passed in {time.perf_counter() - t0:.0f} s; tokens: {tokens['input_tokens']:,} in "
          f"(+{tokens['cache_read_input_tokens']:,} from the cache), {tokens['output_tokens']:,} out")
    if args.out:
        json.dump(rows, open(args.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("answers saved to", args.out)
    return 0 if passed == len(cases) else 1


if __name__ == "__main__":
    sys.exit(main())
