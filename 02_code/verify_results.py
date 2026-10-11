#!/usr/bin/env python3
"""verify_results.py - compares the results.json files of a fresh run (WORK/<step>/results.json) with the reference files in 03_numbers.
Numbers must agree to a relative tolerance of 1e-6 (the pipelines are deterministic); text and lists must match exactly."""
import argparse, json, os, math, sys

PAIRS = [("s2/results.json", "sample2_results.json"), ("s4/results.json", "sample4_wall_results.json"), ("s3/results.json", "sample3_base_results.json"),
         ("s3v/walls.json", "sample3_vessel_walls.json"), ("s1/results.json", "sample1_results.json"), ("s1_pot/results.json", "sample1_pot_results.json")]

def diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b: out.append("%s/%s: missing on one side" % (path, k))
            else: out += diff(a[k], b[k], "%s/%s" % (path, k))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append("%s: list length %d vs %d" % (path, len(a), len(b)))
        else:
            for i, (x, y) in enumerate(zip(a, b)): out += diff(x, y, "%s[%d]" % (path, i))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-9): out.append("%s: %r vs %r" % (path, a, b))
    elif a != b: out.append("%s: %r vs %r" % (path, a, b))
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); ap.add_argument("--reference", required=True); a = ap.parse_args(); bad = 0
    for w, r in PAIRS:
        pw, pr = os.path.join(a.work, w), os.path.join(a.reference, r)
        if not os.path.exists(pw): print("SKIPPED (not run):", w); continue
        d = diff(json.load(open(pw)), json.load(open(pr)))
        print("%-22s %s" % (w, "identical to the reference" if not d else "%d differences" % len(d)))
        for x in d[:8]: print("    ", x)
        bad += len(d)
    sys.exit(1 if bad else 0)
