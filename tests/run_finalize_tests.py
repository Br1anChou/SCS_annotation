#!/usr/bin/env python3
"""Regression tests for finalize_annotation.py (stage-2 finalizer).

Run:  python3 tests/run_finalize_tests.py   (exit 0 = pass, 1 = fail)

Checks that an overrides table:
  * replaces the annotation columns for the clusters it names,
  * recomputes confidence from an overridden score when confidence is omitted,
  * leaves non-overridden clusters and ALL original columns untouched,
  * keeps the partial annotation fields it does not mention (e.g. marker_basis).
"""
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "scripts" / "finalize_annotation.py"
CAND = HERE / "_finalize_candidate.tsv"
OVR = HERE / "_finalize_overrides.tsv"
OUT = HERE / "_finalize_final.tsv"

CANDIDATE = (
    "cluster\tgene\tmean_expression\tcell_ratio\tmanual_type\tscore\tconfidence\t"
    "marker_basis\tevidence_support\tannotation_note\n"
    "0\tPLP1\t4.2\t0.99\tOligodendrocyte\t0.88\thigh\tASPA+CNP\tcand-evidence\tcand-note\n"
    "0\tMBP\t2.6\t0.92\tOligodendrocyte\t0.88\thigh\tASPA+CNP\tcand-evidence\tcand-note\n"
    "1\tAQP4\t2.0\t0.90\tAstrocyte\t0.70\tmoderate\tAQP4\tcand-evidence\tcand-note\n"
)
# Override cluster 0's type + score only (no confidence -> must auto-derive).
OVERRIDES = (
    "cluster\tmanual_type\tscore\n"
    "0\tOligodendrocyte-1\t0.95\n"
)


def main():
    CAND.write_text(CANDIDATE)
    OVR.write_text(OVERRIDES)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(CAND), "--overrides", str(OVR),
         "--output", str(OUT)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        print(proc.stdout, proc.stderr)
        raise SystemExit("finalize script crashed")

    with open(OUT) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        rows = list(reader)
    by_cluster = {}
    for r in rows:
        by_cluster.setdefault(r["cluster"], r)

    checks = []
    def check(name, cond): checks.append((name, bool(cond)))

    c0, c1 = by_cluster["0"], by_cluster["1"]

    check("c0 manual_type overridden", c0["manual_type"] == "Oligodendrocyte-1")
    check("c0 score overridden", c0["score"] == "0.95")
    check("c0 confidence auto-derived from score", c0["confidence"] == "high")
    check("c0 unmentioned field kept (marker_basis)", c0["marker_basis"] == "ASPA+CNP")
    check("c0 original columns preserved", c0["gene"] == "PLP1" and c0["mean_expression"] == "4.2")
    check("c1 untouched (manual_type)", c1["manual_type"] == "Astrocyte")
    check("c1 untouched (score/confidence)", c1["score"] == "0.70" and c1["confidence"] == "moderate")
    check("header preserves original + annotation cols",
          header[:4] == ["cluster", "gene", "mean_expression", "cell_ratio"])
    check("override applied to every row of cluster",
          all(r["manual_type"] == "Oligodendrocyte-1" for r in rows if r["cluster"] == "0"))

    width = max(len(n) for n, _ in checks)
    failed = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{width}}")
        failed += not ok

    for p in (CAND, OVR, OUT):
        p.unlink(missing_ok=True)
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
