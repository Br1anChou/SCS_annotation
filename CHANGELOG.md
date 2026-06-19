# Changelog

All notable changes to the scs-annotator skill are documented here.
This project follows [Semantic Versioning](https://semver.org/).

## [1.0.0]

Major behavioral overhaul of ambient-RNA / doublet handling, plus packaging.

### Changed
- **Data-driven ambient-RNA detection.** `flag_ambient_rna` no longer relies on
  a hard-coded gene list alone; it flags any gene that is broadly present across
  clusters yet never cluster-specific (low specificity everywhere). Catches
  tissue-specific bleeders such as `TAC1` that the seed list missed.
- **Specificity-gated doublet detection.** `detect_doublets` now requires *both*
  conflicting lineages to be cluster-specific (high `cell_ratio` **and** high
  specificity) before calling `Mixed`/`Doublet`. Ambient bleed can no longer
  fake a doublet, while genuine cell-mixing events are still caught.
- **Borderline conflicts preserved.** When one lineage is specific and the other
  is a sub-specific/non-ambient marker, the dominant type is kept but a
  `Doublet QC` note is emitted and confidence is capped — the doublet hypothesis
  is never silently dropped.
- **Label-aware `marker_basis`.** Lists the markers supporting the assigned type
  (panel hits, ambient excluded) instead of the raw top-ranked genes.

### Fixed
- Seed ambient genes (`PLP1`/`MBP`) are no longer flagged when they are the
  strongly-specific defining marker of an oligodendrocyte cluster — which had
  blinded doublet detection in that cluster.

### Added
- `README.md` with full input/output specification.
- `tests/run_tests.py` + `tests/fixture_doublet.tsv` regression test (genuine
  doublet detected, ambient-only bleed rejected, borderline surfaced).
- Cell-level doublet confirmation guidance (Scrublet / DoubletFinder / UMAP) in
  `references/databases.md`; specificity-gated doublet section in
  `references/markers.md`.
- `version` field in `SKILL.md`.

## [0.1.0]

Initial version of the skill.

- Four-dimension candidate scoring: statistical exclusivity, tissue priors,
  doublet/ambient detection, hierarchical (broad class → subtype) annotation.
- Ambient detection limited to a hard-coded seed gene list.
- Doublet detection based on marker `cell_ratio` thresholds (no specificity gate).
- `references/markers.md` marker panels and `references/databases.md` query guide.
