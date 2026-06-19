#!/usr/bin/env python3
"""Regression tests for ambient-aware, specificity-gated doublet detection.

Run:  python3 tests/run_tests.py   (exit 0 = pass, 1 = fail)

The fixture has three clusters that exercise the three behaviours the
ambient/doublet rework must guarantee simultaneously:

  cluster 0  astrocyte with a borderline, NON-ambient neuron marker (SYT1)
             -> dominant Astrocyte call kept, but a soft "Doublet QC" note is
                surfaced (never silently dropped) and confidence is capped.
  cluster 1  genuine oligodendrocyte+neuron doublet: BOTH lineages are
             cluster-specific -> must still be flagged Mixed/Doublet.
  cluster 2  microglia whose only cross-lineage signal is ambient PENK
             -> pure background must NOT raise a doublet flag.
"""
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "score_annotations.py"
FIXTURE = HERE / "fixture_doublet.tsv"
OUT = HERE / "_fixture_doublet.annotated.txt"


def run():
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(FIXTURE),
         "--tissue", "brain", "--output", str(OUT)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr)
        raise SystemExit("script crashed")
    rows_by_cluster = {}
    with open(OUT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows_by_cluster.setdefault(r["cluster"], r)  # one row per cluster
    return proc.stdout, rows_by_cluster


def main():
    stdout, rows = run()
    checks = []

    def check(name, cond):
        checks.append((name, bool(cond)))

    # Ambient detection must catch PENK (broad, low-specificity bleeder).
    check("PENK flagged ambient", "PENK" in stdout and "ambient" in stdout.lower())

    c0, c1, c2 = rows["0"], rows["1"], rows["2"]

    # cluster 1: genuine doublet still detected.
    check("c1 real doublet flagged Mixed/Doublet",
          c1["manual_type"].startswith(("Mixed", "Doublet")))

    # cluster 0: dominant Astrocyte kept, soft conflict surfaced + confidence capped.
    check("c0 not hard-doublet labelled",
          not c0["manual_type"].startswith(("Mixed", "Doublet")))
    check("c0 is astrocyte", "Astro" in c0["manual_type"])
    check("c0 surfaces Doublet QC note", "Doublet QC" in c0["annotation_note"])
    check("c0 confidence not 'high'", c0["confidence"] != "high")

    # cluster 2: pure ambient overlap must NOT raise a doublet flag.
    check("c2 is microglia", "Microglia" in c2["manual_type"])
    check("c2 no false Doublet QC note", "Doublet QC" not in c2["annotation_note"])
    check("c2 clean (not Mixed/Doublet)",
          not c2["manual_type"].startswith(("Mixed", "Doublet")))

    width = max(len(n) for n, _ in checks)
    failed = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{width}}")
        failed += not ok

    OUT.unlink(missing_ok=True)
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
