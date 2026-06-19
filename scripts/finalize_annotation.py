#!/usr/bin/env python3
"""Stage-2 finalizer for scs-annotator.

`score_annotations.py` produces *candidate* annotations. After you validate them
(marker specificity + database/literature review), record your final per-cluster
calls in an overrides table and run this to produce the final annotated file.

It keeps every original column and replaces only the six annotation columns for
the clusters you override, leaving the script's candidate values for the rest —
so a curated, reproducible final table is one command away, with no risk of
mangling the original numeric columns.

Usage:
    python finalize_annotation.py <candidate.tsv> --overrides <overrides.tsv> \
        [--output <final.tsv>]

Inputs:
    candidate  Output of score_annotations.py (original columns + the six
               annotation columns). A raw marker table also works — the six
               annotation columns are appended if absent.
    overrides  TSV / CSV / JSON keyed by `cluster`, holding any subset of the
               six annotation columns. Empty cells / absent columns keep the
               candidate value. If `score` is overridden but `confidence` is
               not, the confidence bucket is recomputed from the score.

The annotation columns are: manual_type, score, confidence, marker_basis,
evidence_support, annotation_note.
"""
import argparse
import csv
import json
from pathlib import Path

ANNO_COLS = ["manual_type", "score", "confidence",
             "marker_basis", "evidence_support", "annotation_note"]


def _delim(path: str) -> str:
    return "," if Path(path).suffix.lower() == ".csv" else "\t"


def _cluster_key(value) -> str:
    """Normalize a cluster id so '0', '0.0', 0 all match."""
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).strip()


def confidence_label(score) -> str:
    """Same thresholds as score_annotations.py (kept inline to stay standalone)."""
    try:
        s = float(score)
    except (TypeError, ValueError):
        return ""
    if s >= 0.80:
        return "high"
    if s >= 0.50:
        return "moderate"
    if s >= 0.30:
        return "low"
    return "unreliable"


def load_overrides(path: str) -> dict[str, dict[str, str]]:
    """cluster_key -> {annotation_col: value} for non-empty provided fields."""
    overrides: dict[str, dict[str, str]] = {}
    if Path(path).suffix.lower() == ".json":
        with open(path) as f:
            raw = json.load(f)
        for key, fields in raw.items():
            picked = {c: str(fields[c]) for c in ANNO_COLS
                      if c in fields and str(fields[c]).strip() != ""}
            if picked:
                overrides.setdefault(_cluster_key(key), {}).update(picked)
        return overrides

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter=_delim(path))
        if "cluster" not in (reader.fieldnames or []):
            raise ValueError("overrides table needs a 'cluster' column")
        for row in reader:
            key = _cluster_key(row.get("cluster", ""))
            if not key:
                continue
            picked = {c: row[c].strip() for c in ANNO_COLS
                      if c in row and row[c] is not None and row[c].strip() != ""}
            if picked:
                overrides.setdefault(key, {}).update(picked)
    return overrides


def main():
    ap = argparse.ArgumentParser(
        description="Apply validated per-cluster annotations onto candidate output")
    ap.add_argument("candidate", help="Candidate annotation file (score_annotations.py output)")
    ap.add_argument("--overrides", required=True,
                    help="TSV/CSV/JSON of final per-cluster calls (keyed by cluster)")
    ap.add_argument("--output", default=None,
                    help="Output path (default: anno_final_annotation.txt next to candidate)")
    args = ap.parse_args()

    out = args.output or str(Path(args.candidate).with_name("anno_final_annotation.txt"))
    if Path(out).resolve() == Path(args.candidate).resolve():
        raise ValueError("Output path resolves to the candidate input; refusing to overwrite")

    overrides = load_overrides(args.overrides)

    with open(args.candidate, newline="") as f:
        reader = csv.DictReader(f, delimiter=_delim(args.candidate))
        header = list(reader.fieldnames or [])
        if "cluster" not in header:
            raise ValueError("candidate file needs a 'cluster' column")
        rows = list(reader)

    out_header = header + [c for c in ANNO_COLS if c not in header]

    candidate_clusters = {_cluster_key(r.get("cluster", "")) for r in rows}
    applied = sorted(k for k in overrides if k in candidate_clusters)
    missing = sorted(k for k in overrides if k not in candidate_clusters)

    with open(out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(out_header)
        for r in rows:
            vals = dict(r)
            for c in ANNO_COLS:
                vals.setdefault(c, "")
            ov = overrides.get(_cluster_key(r.get("cluster", "")), {})
            vals.update(ov)
            # Recompute confidence when the score is overridden but confidence is not.
            if "score" in ov and "confidence" not in ov:
                lab = confidence_label(vals["score"])
                if lab:
                    vals["confidence"] = lab
            writer.writerow([vals.get(c, "") for c in out_header])

    print(f"Wrote {out}")
    print(f"  rows: {len(rows)}  |  clusters overridden: {len(applied)}"
          + (f" ({', '.join(applied)})" if applied else ""))
    if missing:
        print(f"  WARNING: overrides for clusters absent from candidate, ignored: "
              f"{', '.join(missing)}")


if __name__ == "__main__":
    main()
