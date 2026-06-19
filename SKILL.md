---
name: scs-annotator
version: 1.0.0
description: >-
  Rigorous single-cell RNA-seq cluster cell-type annotation. Use when given a cluster marker table or Seurat FindAllMarkers output (columns like cluster/gene/avg_log2FC/pct.1 or cluster/gene/mean_expression/cell_ratio) and asked to assign cell types, name clusters, validate markers, or produce publication-style annotations. Distinguishes ambient-RNA bleed from genuine doublets via marker specificity, and flags candidate doublets/mixed clusters for per-cell QC. Triggers on: 单细胞注释, 细胞注释, 细胞类型鉴定, 细胞类型注释, marker基因注释, 注释marker表, 双细胞/doublet 判定, annotate clusters, cell type identification.
---

# scs-annotator

Annotate marker tables by preserving original columns and appending:

`manual_type`, `score`, `confidence`, `marker_basis`, `evidence_support`, `annotation_note`

Use the exact column name `manual_type` (not `mannual_type`). Default output for downstream analysis is `anno_candidate_annotation.txt` in the input directory. Default species is human; use `--species mouse` only when explicitly told. The bundled script is a candidate generator, not final literature validation. Publication-facing results require database and paper-level evidence.

## Input

Preferred TSV: `cluster | gene | cell_type | mean_expression | cell_count | cell_ratio`

Also accepts Seurat-like tables with `cluster`, `gene`, `avg_log2FC` or `avg_logFC`, and `pct.1`.

## Naming

- One unique, concise `manual_type` per cluster; repeat it on every marker row in that cluster.
- Use classic names: `Oligodendrocyte`, `Astrocyte`, `OPC`, `Microglia`, `T-cell`, `B-cell`, `Endothelial`, `D1-MSN`, `D2-MSN`, `Interneuron`.
- Add `(geneA+geneB)` only to distinguish clusters that would otherwise have the same name: `D1-MSN(PDYN+TAC1)`, `Oligodendrocyte(ASPA+CNP)`.
- If gene suffixes still collide, add `-1/-2/-3` after the base name: `Oligodendrocyte-1(ASPA+CNP)`, `Oligodendrocyte-2(ASPA+CNP)`.
- Keep generic lineage genes out of suffixes when possible: `SNAP25`, `SYT1`, `MBP`, `PLP1`, `PTPRC`, `PECAM1`.

## Decision Rules

1. **Statistical exclusivity first.** Prefer markers with `cell_ratio > 0.60` and cluster-specific `mean_expression`; downweight genes with similar expression across many clusters.
2. **Use tissue priors.** Striatum/subcortical: `TAC1/DRD1/PDYN` for D1-MSN and `PENK/DRD2/GPR6/ADORA2A` for D2-MSN. Cortex: `SLC17A7/SATB2/CAMK2A` for excitatory neurons and `GAD1/GAD2/SST/PVALB/VIP` for inhibitory neurons. Multi-region brain keeps both priors active.
3. **Be conservative with neurons.** Call neuron subtypes only when subtype markers dominate. Generic neuronal markers alone (`SNAP25`, `SYT1`, etc.) are not enough for D1/D2/interneuron subtype calls.
4. **Prioritize glial lineage markers.** If lineage glial markers dominate (`SOX10/OLIG1/OLIG2` for oligodendrocyte lineage; `PDGFRA/GPR17` for OPC; `AQP4/SLC1A3/GFAP` for astrocyte; `TMEM119/P2RY12/CX3CR1` for microglia), prefer glia even if weak neuron markers appear.
5. **Mixed and artifacts — separate ambient bleed from real doublets by specificity.** Call `Mixed(Glia+Neuron)` / `Doublet` only when **both** conflicting lineages are *cluster-specific* (high `cell_ratio` **and** high specificity vs other clusters), e.g. `PLP1/MBP + SNAP25/SYT1` both enriched here and low elsewhere. A second lineage that is high-ratio but globally present / low-specificity (e.g. striatal `PENK`, myelin `PLP1/MBP` bleeding everywhere) is **ambient RNA, not a doublet** — assign the dominant type. But never silently erase the possibility: when any cross-lineage signal remains, keep it in `annotation_note`, lower `confidence`, and recommend per-cell doublet confirmation. Cluster-level marker conflict is only a *candidate* doublet flag.
6. **Layer annotation.** Broad class first: neuron, glia, immune, vascular, progenitor. Then subtype using specific marker panels and tissue context.

## Evidence Workflow

- Actively web-search ambiguous, low-confidence, neuron-subtype, mixed/doublet, rare, or publication-facing annotations.
- **Confirm candidate doublets at the cell level.** A cluster-level marker conflict cannot prove a doublet — verify with per-cell doublet scores (Scrublet / DoubletFinder) and UMAP position (doublet clusters often sit between two parents) before committing to a `Doublet`/`Mixed` call, and before dismissing one.
- Use at least two sources when possible: PanglaoDB, CellMarker 2.0, Allen Brain Atlas/Cell Types, and PubMed primary literature.
- Do not claim paper-level evidence from the script alone. If evidence is missing or conflicting, lower `score/confidence` and say so in `annotation_note`.
- For final manuscripts, every cluster-level `manual_type` should be backed by marker specificity plus database/literature evidence.

## Resources

- `scripts/score_annotations.py`: deterministic candidate annotation and scoring. Flags ambient RNA data-drivenly (broadly-present, never-specific genes) and gates doublet calls by specificity; emits soft `Doublet QC` notes for borderline conflicts.
- `references/markers.md`: compact marker panel.
- `references/databases.md`: search/query guide, including per-cell doublet confirmation.
- `tests/run_tests.py` + `tests/fixture_doublet.tsv`: regression test locking in real-doublet detection vs ambient-bleed rejection. Run `python3 tests/run_tests.py` after editing the scoring script.
