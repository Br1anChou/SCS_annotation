# scs-annotator-skill

Rigorous single-cell RNA-seq cluster cell-type annotation, packaged as a
[Claude Code](https://claude.com/claude-code) skill. It takes a per-cluster
marker table, **preserves every original column**, and appends six annotation
columns describing the called cell type and the evidence behind it.

The bundled script is a **candidate generator** — a deterministic, reproducible
first pass. It does *not* replace database/literature validation; publication
calls still require that step (see [Caveats](#caveats)).

---

## Contents

- [Install](#install-as-a-claude-code-skill)
- [Input](#input) — formats, columns, examples
- [Running the script](#running-the-script) — CLI, flags, stdout
- [Output](#output) — appended columns, value ranges, examples
- [Finalizing (stage 2)](#finalizing-stage-2) — apply your validated calls
- [How calls are made](#how-calls-are-made) — ambient RNA & doublets
- [Tests](#tests) · [Caveats](#caveats) · [Layout](#repository-layout)

---

## Install (as a Claude Code skill)

```bash
git clone https://github.com/Br1anChou/SCS_annotation.git ~/.claude/skills/scs-annotator
```

Invoke in Claude Code with `/scs-annotator`, or let it trigger on phrasings like
"单细胞注释 / annotate clusters / cell type identification". You can also run the
script standalone (no Claude Code required) — see [Running the script](#running-the-script).

---

## Input

A **tab-separated** (`.tsv`/`.txt`) or **comma-separated** (`.csv`) table with a
header row. The delimiter is auto-detected from the file extension (`.csv` → comma,
everything else → tab). Each row is **one marker gene in one cluster**, so a
cluster appears across many rows.

Two header layouts are accepted.

### Format A — normalized expression table (preferred)

| Column            | Type   | Required | Description |
|-------------------|--------|----------|-------------|
| `cluster`         | int    | **yes**  | Cluster ID. Repeats across the gene rows belonging to it. |
| `gene`            | string | **yes**  | Gene symbol. Case-insensitive; matched upper-cased against the panels (`Aqp4` = `AQP4`). |
| `cell_type`       | string | no       | Optional prior label (e.g. from a reference tool). Used only as a **weak hint**; the call is driven by the markers, not this column. Preserved in the output. |
| `mean_expression` | float  | yes\*    | Mean expression of the gene in the cluster. Used for specificity (fold-change vs. other clusters). |
| `cell_count`      | int    | yes\*    | Number of cells in the cluster expressing the gene. Carried through; not central to scoring. |
| `cell_ratio`      | float  | yes\*    | Fraction of cells in the cluster expressing the gene (0–1). The primary statistical signal. |

\* Format A requires **all three** of `mean_expression`, `cell_count`, `cell_ratio`.

**Example (`.tsv`):**

```
cluster	gene	cell_type	mean_expression	cell_count	cell_ratio
0	PLP1	Oligodendrocyte	4.2035	23535	0.9969
0	MBP	Oligodendrocyte	2.6276	21732	0.9205
0	ASPA	Oligodendrocyte	2.2191	18976	0.8038
0	CNP	Oligodendrocyte	1.4070	15132	0.6409
1	DRD1	D1_MSN	0.8909	11064	0.7384
1	PDYN	D1_MSN	0.6722	8699	0.5806
```

### Format B — Seurat `FindAllMarkers` output

If the normalized columns are absent, the script accepts a Seurat marker table
that has `cluster`, `gene`, **`pct.1`**, and one of **`avg_log2FC`** / **`avg_logFC`**.

| Seurat column            | Mapped to         | Notes |
|--------------------------|-------------------|-------|
| `cluster`                | `cluster`         | required |
| `gene`                   | `gene`            | required |
| `avg_log2FC`/`avg_logFC` | `mean_expression` | negative values clamped to 0 |
| `pct.1`                  | `cell_ratio`      | fraction expressing in-cluster |
| `cell_count` (if present)| `cell_count`      | else 0 |

```
cluster	gene	avg_log2FC	pct.1	pct.2	p_val_adj
0	PLP1	3.91	0.997	0.21	0.0
0	MBP	2.88	0.920	0.18	0.0
1	DRD1	2.10	0.738	0.05	1e-300
```

> **Tip — specificity needs all clusters.** Specificity = a gene's expression in
> its cluster ÷ its highest expression in any *other* cluster. So include every
> cluster's markers in one file (export the full `FindAllMarkers` table rather
> than a single cluster) — otherwise ambient/doublet detection has nothing to
> compare against.

---

## Running the script

```bash
python3 scripts/score_annotations.py <input.tsv> [options]
```

| Flag             | Default                         | Description |
|------------------|---------------------------------|-------------|
| `<input>`        | —                               | Path to the marker table (positional, required). |
| `--tissue`       | auto-detected from markers      | Tissue context, e.g. `brain`, `striatum`, `cortex`, `immune`, `liver`. Drives tissue priors and the brain glia/neuron rules. |
| `--species`      | `human`                         | `human` or `mouse` (symbols are upper-cased for matching either way). |
| `--output`       | `anno_candidate_annotation.txt` next to the input | Output path. **Refuses to overwrite the input file.** |
| `--markers`      | —                               | Extra marker panel as TSV (`cell_type<TAB>GENE1<TAB>GENE2…`). |
| `--markers-json` | —                               | Extra marker panel as JSON (`{"cell_type": ["GENE1", …]}`). |

**Example:**

```bash
python3 scripts/score_annotations.py markers.tsv --tissue brain --species human
```

**What it prints to stdout** (diagnostics; the annotated table goes to the output file):

```
# Inferred tissue: brain
# Flagged ambient RNA genes: MBP, PENK, PLP1, TAC1

## Annotation Summary (13 clusters, tissue: brain, species: human)
Cl   Original          Broad     manual_type                Score  Conf      Top Markers
----------------------------------------------------------------------------------------
0    Oligodendrocyte   Glia      Oligodendrocyte(ASPA+CNP)  0.88   high      MBP, PLP1, ASPA
1    D1_MSN            Neuron    D1-MSN(DRD1+RELN)          0.73   moderate  RELN, DRD1, PDYN
...
```

- `# Inferred tissue` — the tissue prior actually used.
- `# Flagged ambient RNA genes` — genes treated as background bleed (excluded
  from doublet/lineage decisions). See [How calls are made](#how-calls-are-made).

---

## Output

A table written to `--output` (default `anno_candidate_annotation.txt` beside the
input). It is **all your original columns, unchanged**, with **six columns appended**.
(Any pre-existing annotation columns from a previous run — including the legacy
misspelling `mannual_type` — are stripped and regenerated, so re-runs stay clean.)

The same cluster-level call is **repeated on every marker row of that cluster**,
so you can group by `cluster` or filter rows freely without losing the label.

### Appended columns

| Column             | Type   | Range / values | Meaning |
|--------------------|--------|----------------|---------|
| `manual_type`      | string | e.g. `Astrocyte`, `D1-MSN(DRD1+RELN)`, `Mixed(Glia+Neuron)` | The called cell type — one concise, unique label per cluster. A `(geneA+geneB)` suffix is added only to disambiguate clusters that would otherwise share a name; `-1/-2` ordinals if suffixes still collide. |
| `score`            | float  | `0.00`–`1.00`  | Heuristic confidence. Sum of: statistical support (≤0.30) + canonical-marker match (≤0.30) + lineage purity (≤0.25) + subtype marker strength (≤0.15). Database/literature concordance is **not** awarded by the script, so a perfectly real type still tops out below 1.0. |
| `confidence`       | string | `high` ≥0.80 · `moderate` ≥0.50 · `low` ≥0.30 · `unreliable` <0.30 | Bucketed `score`. |
| `marker_basis`     | string | `+`-joined genes, e.g. `ASPA+CNP` | The markers actually supporting the call (panel hits in rank order, ambient excluded) — not just the top-ranked genes. |
| `evidence_support` | string | free text | What the call rests on and where to verify it (named databases/atlases). Always phrased as a **candidate** needing validation. |
| `annotation_note`  | string | free text | Caveats and reasoning: ambient/background warnings, glia-vs-neuron dominance, and — for borderline clusters — a `\| Doublet QC: …` suffix telling you to confirm per-cell. |

### Example output rows

A clean, high-confidence glial cluster (original columns 1–6 preserved, 7–12 appended):

```
cluster  gene  cell_type        mean_expression  cell_count  cell_ratio  →  manual_type                score  confidence  marker_basis  evidence_support                                                              annotation_note
0        PLP1  Oligodendrocyte  4.2035           23535       0.9969      →  Oligodendrocyte(ASPA+CNP)  0.88   high        ASPA+CNP      Candidate: MBP/PLP1/ASPA/CNP glial marker panel; validate with CellMarker 2.0, PanglaoDB, Allen/PubMed  Glial lineage markers dominate; weaker neuronal markers should be treated as possible background
```

A genuine doublet (both lineages cluster-specific) and a borderline cluster
flagged for QC:

```
manual_type         score  confidence  annotation_note
Mixed(Glia+Neuron)  0.30   low         Low score; use as candidate annotation pending manual/literature review
Astrocyte           0.65   moderate    Glial lineage markers dominate … | Doublet QC: cross-lineage markers co-occur (Glia+Neuron) but the secondary lineage is sub-specific/ambient — bleed most likely; rule out doublet with per-cell Scrublet/DoubletFinder + UMAP position
```

Read it as: `manual_type` is the call; `score`/`confidence` is how much to trust
it; `marker_basis` is *why*; `evidence_support` tells you which database to check;
`annotation_note` is the fine print. **A `Doublet QC:` note means the dominant
type was assigned but a doublet is not excluded — confirm at the cell level.**

---

## Finalizing (stage 2)

`score_annotations.py` produces **candidates**. The intended workflow is two
stages:

```
markers.tsv ──score_annotations.py──▶ candidate ──(you validate)──▶ overrides.tsv
                                          │                              │
                                          └────finalize_annotation.py────┴──▶ final
```

After validating the candidates (marker specificity + database/literature),
record your final per-cluster calls in an **overrides** table and apply them with
`scripts/finalize_annotation.py`. It keeps every original column and replaces only
the annotation columns for the clusters you name, leaving the script's candidate
values for the rest — so you never hand-edit the numeric columns.

```bash
python3 scripts/finalize_annotation.py anno_candidate_annotation.txt \
    --overrides my_overrides.tsv --output anno_final_annotation.txt
```

**Overrides table** — TSV/CSV/JSON keyed by `cluster`, with any subset of the six
annotation columns. Empty cells / absent columns keep the candidate value. If you
override `score` but omit `confidence`, the confidence bucket is recomputed for you.

```
cluster	manual_type	score	annotation_note
3	Astrocyte-1	0.90	Clean astrocyte; weak neuron signal is ambient PENK/TAC1, not a doublet.
8	Interneuron	0.70	GABAergic interneuron: LHX6/SST/NPY/VIP enriched; candidate D1-MSN overridden.
```

JSON form: `{"3": {"manual_type": "Astrocyte-1", "score": "0.90"}, ...}`.

A complete worked example for the striatal dataset is in
[`examples/striatum_overrides.tsv`](examples/striatum_overrides.tsv) (all 13
clusters, with marker basis, evidence, and notes).

---

## How calls are made

Full decision rules live in [`SKILL.md`](SKILL.md); the two that matter most for
reading the output:

**Ambient RNA vs. doublet (specificity-gated).** A gene that is highly expressed
but spread across most clusters at low specificity (e.g. striatal `PENK`, myelin
`PLP1`/`MBP`) is treated as **ambient bleed** and excluded from lineage/doublet
decisions — it will not, by itself, turn a clean cluster into a doublet. A
`Mixed`/`Doublet` call is made only when **both** conflicting lineages are
*cluster-specific* (high `cell_ratio` **and** high specificity). Borderline cases
are never silently dropped: the dominant type is kept, confidence is capped, and a
`Doublet QC` note is emitted. See [`references/databases.md`](references/databases.md)
for per-cell confirmation (Scrublet / DoubletFinder / UMAP).

**Glia-first in brain tissue.** When definitive glial lineage markers dominate and
neuronal markers are weak, the cluster is called glia even if faint neuron markers
appear (they are usually background).

---

## Tests

```bash
python3 tests/run_tests.py           # scorer: ambient vs doublet   (exit 0 = pass)
python3 tests/run_finalize_tests.py  # finalizer: overrides merge   (exit 0 = pass)
```

`run_tests.py` ([`tests/fixture_doublet.tsv`](tests/fixture_doublet.tsv)) locks in
three scorer behaviours at once: a genuine doublet is detected, ambient-only bleed
is **not** flagged as a doublet, and a borderline secondary lineage is surfaced
(not dropped). `run_finalize_tests.py` checks that overrides replace only the named
clusters/fields, recompute confidence from an overridden score, and preserve every
original column. Run both after editing the scripts.

---

## Caveats

- The script is a **candidate generator**, not final validation. It does not award
  database/literature concordance — every `manual_type` should be cross-checked
  against PanglaoDB, CellMarker 2.0, Allen Brain Atlas, and primary literature
  before publication.
- Cluster-level marker conflict cannot prove a doublet; confirm at the cell level.
- The optional `cell_type` input column is a *weak hint only* — a wrong prior label
  will not override clear marker evidence, but don't rely on it being authoritative.

---

## Repository layout

```
SKILL.md                          # skill definition + decision rules (entry point)
scripts/score_annotations.py      # stage 1: deterministic candidate annotation & scoring
scripts/finalize_annotation.py    # stage 2: apply validated per-cluster calls onto candidates
references/markers.md             # canonical marker panels by tissue / cell type
references/databases.md           # marker-DB query guide + doublet confirmation
examples/striatum_overrides.tsv   # worked overrides example (striatal dataset)
tests/run_tests.py                # regression test (ambient vs doublet)
tests/run_finalize_tests.py       # regression test (finalizer)
tests/fixture_doublet.tsv         # test fixture
CHANGELOG.md                      # version history
```
