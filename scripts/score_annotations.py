#!/usr/bin/env python3
"""
scs-annotator: Four-dimension single-cell annotation scoring.

Dimension 1: Statistical exclusivity (cell_ratio, mean_expression specificity)
Dimension 2: Tissue-anatomical priors (tissue-specific marker weighting)
Dimension 3: Doublet & ambient RNA detection
Dimension 4: Hierarchical annotation (broad class → subtype)

Usage:
    python score_annotations.py <input.tsv> --tissue brain --species human
"""

import json
import argparse
import re
import csv
from collections import defaultdict
from pathlib import Path


# ── Broad class definitions (Dimension 4) ──

BROAD_CLASSES = {
    "Neuron": {
        "generic": {"SNAP25", "SYT1", "RBFOX3", "MAP2", "STMN2", "ELAVL4", "NSG1", "NEFL"},
        "excitatory": {"SLC17A7", "SLC17A6", "CAMK2A", "SATB2", "BCL11B"},
        "inhibitory": {"GAD1", "GAD2", "SLC32A1"},
        "msn": {"PPP1R1B", "DRD1", "DRD2", "PDYN", "TAC1", "PENK", "GPR6", "ADORA2A",
                "GPR149", "OTOF", "CPNE4", "SPHKAP", "RELN", "ISL1", "OPRM1"},
        "interneuron": {"PVALB", "SST", "VIP", "NPY", "CCK", "CALB2", "HTR3A",
                        "LHX6", "KCNS3", "KCNC1", "KCNC2", "CHAT"},
        "cortical": {"CUX1", "CUX2", "FEZF2", "RORB", "TLE4", "THEMIS"},
    },
    "Glia": {
        "oligodendrocyte": {"MBP", "PLP1", "MOG", "SOX10", "CLDN11", "CNP", "ASPA"},
        "opc": {"PDGFRA", "OLIG1", "OLIG2", "CSPG4", "GPR17"},
        "astrocyte": {"AQP4", "GFAP", "SLC1A3", "ALDH1L1", "GJA1"},
        "microglia": {"TMEM119", "P2RY12", "CX3CR1", "CSF1R", "ITGAM"},
    },
    "Immune": {
        "t_cell": {"CD3E", "CD3D", "CD3G", "CD2", "CD8A", "CD4"},
        "b_cell": {"CD19", "MS4A1", "CD79A", "CD79B", "PAX5"},
        "nk_cell": {"NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"},
        "myeloid": {"CD14", "CD68", "CD163", "CSF1R", "ITGAM", "S100A8", "LYZ"},
    },
    "Vascular": {
        "endothelial": {"PECAM1", "CDH5", "VWF", "ESAM", "CLDN5", "ENG"},
        "pericyte": {"PDGFRB", "RGS5", "CSPG4", "ANPEP"},
        "vsmc": {"ACTA2", "MYH11", "TAGLN", "CNN1"},
    },
    "Progenitor": {
        "proliferating": {"MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "BIRC5"},
        "neuroblast": {"DCX", "SOX4", "SOX11", "STMN1"},
        "npc": {"SOX2", "NES", "PAX6", "HES5"},
    },
}

# Marker panels used for conservative brain annotation and QC flags.
BRAIN_TYPE_MARKERS = {
    "Oligodendrocyte": {"PLP1", "MBP", "ASPA", "CNP", "MOG", "MAG", "CLDN11", "MOBP", "OPALIN", "SOX10"},
    "OPC": {"PDGFRA", "GPR17", "OLIG1", "OLIG2", "CSPG4", "SOX10"},
    "Astrocyte": {"AQP4", "SLC1A3", "GFAP", "ALDH1L1", "GJA1", "GPC5", "S100B", "GLUL"},
    "Microglia": {"TMEM119", "P2RY12", "CX3CR1", "CSF1R", "AIF1", "CD68", "ITGAM"},
    "D1-MSN": {"DRD1", "TAC1", "PDYN", "SPHKAP", "CPNE4", "RELN", "GPR149", "OTOF", "PPP1R1B"},
    "D2-MSN": {"PENK", "DRD2", "ADORA2A", "GPR6", "GPR52", "PPP1R1B"},
    "Interneuron": {"GAD1", "GAD2", "SLC32A1", "SST", "NPY", "PVALB", "VIP", "LHX6", "CCK", "CALB2"},
    "Excitatory-neuron": {"SLC17A7", "SLC17A6", "CAMK2A", "SATB2", "CUX1", "CUX2", "BCL11B"},
    "Neuroblast": {"DCX", "SOX4", "SOX11", "STMN1", "ELAVL4", "TUBB3"},
    "Endothelial cell": {"ESAM", "VWF", "PECAM1", "CDH5", "CLDN5", "ENG"},
    "Pericyte": {"PDGFRB", "RGS5", "CSPG4", "ANPEP"},
    "Vascular smooth muscle cell": {"ACTA2", "MYH11", "TAGLN", "CNN1"},
    "T-cell": {"CD3E", "CD3D", "CD3G", "CD2", "CD4", "CD8A"},
    "B-cell": {"CD19", "MS4A1", "CD79A", "CD79B", "PAX5"},
    "NK-cell": {"NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"},
}

GLIA_CONFLICT_MARKERS = (
    BRAIN_TYPE_MARKERS["Oligodendrocyte"] |
    BRAIN_TYPE_MARKERS["OPC"] |
    BRAIN_TYPE_MARKERS["Astrocyte"] |
    BRAIN_TYPE_MARKERS["Microglia"]
)

NEURON_CONFLICT_MARKERS = (
    BROAD_CLASSES["Neuron"]["generic"] |
    BROAD_CLASSES["Neuron"]["excitatory"] |
    BROAD_CLASSES["Neuron"]["inhibitory"] |
    BRAIN_TYPE_MARKERS["D1-MSN"] |
    BRAIN_TYPE_MARKERS["D2-MSN"] |
    BRAIN_TYPE_MARKERS["Interneuron"] |
    BRAIN_TYPE_MARKERS["Neuroblast"]
)

# Mutually exclusive gene pairs: high co-signal supports Mixed/Doublet.
MUTUALLY_EXCLUSIVE = [
    # Glia vs Neuron
    (GLIA_CONFLICT_MARKERS, NEURON_CONFLICT_MARKERS),
    # T-cell vs B-cell
    ({"CD3E", "CD3D", "CD3G"}, {"CD19", "MS4A1", "CD79A", "CD79B"}),
    # Endothelial vs Neuron
    ({"PECAM1", "CDH5", "VWF"}, {"SNAP25", "SYT1", "SLC17A7", "GAD1"}),
    # Fibroblast vs parenchymal
    ({"COL1A1", "COL1A2", "DCN", "LUM"}, {"SNAP25", "SYT1", "MBP", "PLP1"}),
    # Smooth muscle vs others
    ({"ACTA2", "MYH11", "TAGLN"}, {"SNAP25", "MBP", "PECAM1"}),
]

# High-abundance genes prone to ambient RNA contamination (brain)
AMBIENT_PRONE = {
    "PENK",   # Striatum — extremely high expression, contaminates all droplets
    "MBP",    # White matter-rich samples
    "PLP1",   # Very high in oligodendrocytes, bleeds everywhere
    "MALAT1", # Ubiquitously high lncRNA
    "NEAT1",  # Ubiquitously high lncRNA
    "XIST",   # Female samples
}

CLASSIC_NAME_ALIASES = {
    "glia": "Glia",
    "neuron": "Neuron",
    "D1_MSN": "D1-MSN",
    "D2_MSN": "D2-MSN",
    "T-cell": "T-cell",
    "CD4+ T-cell": "CD4+ T-cell",
    "CD8+ T-cell": "CD8+ T-cell",
    "B-cell": "B-cell",
    "Plasma-cell": "Plasma cell",
    "NK-cell": "NK-cell",
    "Dendritic-cell": "Dendritic cell",
    "Endothelial": "Endothelial cell",
    "VSMC": "Vascular smooth muscle cell",
    "NPC": "Neural progenitor cell",
    "NB": "Neuroblast",
    "Interneuron": "Interneuron",
}

CANDIDATE_LABEL_BROAD_CLASS = {
    "D1-MSN": "Neuron",
    "D2-MSN": "Neuron",
    "Interneuron": "Neuron",
    "SST-interneuron": "Neuron",
    "VIP-interneuron": "Neuron",
    "PV-interneuron": "Neuron",
    "Excitatory-neuron": "Neuron",
    "Inhibitory-neuron": "Neuron",
    "Oligodendrocyte": "Glia",
    "Astrocyte": "Glia",
    "OPC": "Glia",
    "Microglia": "Glia",
    "Endothelial cell": "Vascular",
    "Neuroblast": "Progenitor",
}

LOW_VALUE_QUALIFIER_MARKERS = {
    "SNAP25", "SYT1", "RBFOX3", "MAP2", "STMN2", "ELAVL4", "NSG1", "NEFL",
    "MBP", "PLP1", "PTPRC", "PECAM1", "CDH5", "VWF", "MALAT1", "NEAT1",
}

STATE_SIGNATURES = [
    ("proliferating", {"MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "BIRC5"}),
    ("homeostatic", {"TMEM119", "P2RY12", "CX3CR1"}),
    ("activated", {"APOE", "SPP1", "LPL", "AIF1", "CD68"}),
    ("myelinating", {"MOG", "MAG", "MOBP", "OPALIN", "CLDN11"}),
    ("interferon-high", {"ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1"}),
    ("cytotoxic", {"NKG7", "GNLY", "GZMB", "PRF1"}),
]


def _as_float(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _as_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def _first_value(row: dict[str, str], names: list[str], default: str = "") -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return default


def _normalize_gene(gene: str, species: str) -> str:
    """Normalize gene symbols for matching against bundled uppercase markers."""
    gene = gene.strip()
    if species == "mouse":
        return gene.upper()
    return gene.upper()


def _detect_delimiter(path: str) -> str:
    suffix = Path(path).suffix.lower()
    return "," if suffix == ".csv" else "\t"


def load_input(path: str, species: str = "human") -> tuple[list[dict], list[str]]:
    rows = []
    delimiter = _detect_delimiter(path)
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("Input file has no header row")
        header = reader.fieldnames
        required = {"cluster", "gene"}
        missing = required - set(header)
        if missing:
            raise ValueError(
                f"Input must contain at least cluster and gene columns; missing: {sorted(missing)}"
            )

        has_normalized = {"mean_expression", "cell_count", "cell_ratio"}.issubset(header)
        has_seurat = any(c in header for c in ("avg_log2FC", "avg_logFC")) and "pct.1" in header
        if not has_normalized and not has_seurat:
            raise ValueError(
                "Unsupported marker table. Provide normalized columns "
                "(mean_expression, cell_count, cell_ratio) or Seurat columns "
                "(avg_log2FC/avg_logFC and pct.1)."
            )

        for raw in reader:
            if not raw.get("cluster") or not raw.get("gene"):
                continue
            gene = _normalize_gene(raw["gene"], species)
            if has_normalized:
                mean_expression = _as_float(raw.get("mean_expression"))
                cell_count = _as_int(raw.get("cell_count"))
                cell_ratio = _as_float(raw.get("cell_ratio"))
                cell_type = raw.get("cell_type", "unknown") or "unknown"
            else:
                avg_fc = _as_float(_first_value(raw, ["avg_log2FC", "avg_logFC"]))
                mean_expression = max(avg_fc, 0.0)
                cell_count = _as_int(raw.get("cell_count"))
                cell_ratio = _as_float(raw.get("pct.1"))
                cell_type = raw.get("cell_type", "FindAllMarkers") or "FindAllMarkers"

            rows.append({
                "cluster": int(float(raw["cluster"])),
                "gene": gene,
                "cell_type": cell_type,
                "mean_expression": mean_expression,
                "cell_count": cell_count,
                "cell_ratio": cell_ratio,
                "_raw": raw,
            })
    return rows, header


def load_markers(path: str) -> dict[str, set[str]]:
    markers = defaultdict(set)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                for gene in parts[1:]:
                    gene = gene.strip()
                    if gene:
                        markers[parts[0].strip()].add(gene)
    return dict(markers)


def parse_markers_md(path: str) -> dict[str, set[str]]:
    """Parse bundled markers.md into cell_type -> marker set."""
    markers = defaultdict(set)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("| ") and "|" in line[2:]:
                parts = [p.strip() for p in line.split("|")[1:-1]]
                if not parts or parts[0] in ("Cell Type", "Category", "---", "----"):
                    continue
                ct = parts[0]
                genes_str = " ".join(parts[1:]) if len(parts) > 1 else ""
                genes = set(re.findall(r'[A-Z][A-Z0-9]+', genes_str))
                if genes:
                    markers[ct].update(genes)
    return dict(markers)


# ═══════════════════════════════════════════════════════════════
# Dimension 1: Statistical Exclusivity Filter
# ═══════════════════════════════════════════════════════════════

def compute_cluster_stats(rows: list[dict]) -> dict[int, dict]:
    """
    For each cluster, compute per-gene specificity and rank genes
    by cell_ratio × specificity_score.
    """
    # Gather per-cluster per-gene mean_expression
    cluster_expr = defaultdict(dict)
    for row in rows:
        cluster_expr[row["cluster"]][row["gene"]] = row["mean_expression"]

    clusters = sorted(cluster_expr.keys())

    stats = {}
    for cid in clusters:
        gene_info = {}
        for row in rows:
            if row["cluster"] != cid:
                continue
            gene = row["gene"]
            expr_here = row["mean_expression"]
            # Compute expression in other clusters
            other_exprs = [
                cluster_expr[oc].get(gene, 0.0)
                for oc in clusters if oc != cid
            ]
            max_other = max(other_exprs) if other_exprs else 0.001
            # Specificity: fold-change vs next highest cluster
            specificity = expr_here / max_other if max_other > 0 else 10.0
            # rank_score = cell_ratio × specificity
            rank_score = row["cell_ratio"] * min(specificity, 10.0)  # cap at 10x

            gene_info[gene] = {
                "cell_ratio": row["cell_ratio"],
                "mean_expression": expr_here,
                "specificity": round(specificity, 2),
                "passes_ratio_filter": row["cell_ratio"] > 0.60,
                "rank_score": round(rank_score, 4),
            }

        # Sort genes by rank_score descending
        sorted_genes = sorted(gene_info.items(), key=lambda x: -x[1]["rank_score"])
        stats[cid] = {
            "genes_ranked": [g for g, _ in sorted_genes],
            "gene_info": gene_info,
            "top_specific_genes": [g for g, _ in sorted_genes[:10]],
        }
    return stats


def dimension1_score(cluster_id: int, stats: dict, gene: str) -> float:
    """Score 0-1 for how statistically valid a gene is as a marker."""
    info = stats[cluster_id]["gene_info"].get(gene, {})
    if not info:
        return 0.0
    # Must pass ratio filter
    if not info["passes_ratio_filter"]:
        return info["cell_ratio"] * 0.5  # partial credit below 0.60
    specificity = info["specificity"]
    if specificity >= 5:
        return 1.0
    elif specificity >= 3:
        return 0.8
    elif specificity >= 2:
        return 0.5
    else:
        return 0.2


# ═══════════════════════════════════════════════════════════════
# Dimension 2: Tissue Priors
# ═══════════════════════════════════════════════════════════════

TISSUE_SIGNATURES = {
    "striatum": {"DRD1", "DRD2", "TAC1", "PENK", "PDYN", "GPR6", "ADORA2A", "PPP1R1B"},
    "cortex": {"SLC17A7", "SATB2", "BCL11B", "GAD1", "PVALB", "SST", "VIP", "CUX2"},
    "hippocampus": {"PROX1", "CBLN3", "SATB2", "SLC17A6", "ZBTB20"},
    "cerebellum": {"PCP2", "CALB1", "GABRA6", "RGS8", "ITPR1"},
    "midbrain": {"TH", "SLC6A3", "NR4A2", "PITX3", "EN1"},
    "spinal_cord": {"MNX1", "CHAT", "SLC5A7", "NTRK1", "TAC1"},
    "liver": {"ALB", "HNF4A", "ASGR1", "CYP3A4", "KRT19"},
    "kidney": {"NPHS1", "NPHS2", "SLC34A1", "SLC12A1", "UMOD"},
    "lung": {"SFTPC", "SFTPB", "AGER", "SCGB1A1", "FOXJ1"},
    "pancreas": {"GCG", "INS", "SST", "CPA1", "KRT19"},
    "heart": {"MYL7", "MYL2", "MYH7", "MYH6", "NPPA"},
    "immune": {"PTPRC", "CD3E", "CD19", "CD14", "NKG7", "MS4A1"},
}


def infer_tissue(rows: list[dict]) -> str:
    """Infer tissue context from marker genes."""
    all_genes = {r["gene"] for r in rows}
    scores = {}
    for tissue, markers in TISSUE_SIGNATURES.items():
        hits = all_genes & markers
        scores[tissue] = len(hits) / len(markers)
    best = max(scores, key=scores.get)
    if scores[best] < 0.05:
        return "unknown"
    return best


# ═══════════════════════════════════════════════════════════════
# Dimension 3: Doublet & Ambient RNA Detection
# ═══════════════════════════════════════════════════════════════

def detect_doublets(
    cluster_genes: dict[int, list[str]],
    cluster_data: dict[int, list[dict]],
    ambient_genes: set[str] | None = None,
) -> dict[int, str]:
    """
    Check each cluster for mutually exclusive marker co-expression.
    Returns {cluster_id: 'clean' | 'Mixed(lineageA+lineageB)' | 'Doublet'}.
    """
    ambient_genes = ambient_genes or set()
    results = {}
    for cid, genes_list in cluster_genes.items():
        # Get gene -> cell_ratio mapping
        gene_ratio = {}
        for row in cluster_data.get(cid, []):
            gene_ratio[row["gene"]] = row["cell_ratio"]

        top_genes = set(genes_list[:15]) - ambient_genes
        conflicts = []
        for set_a, set_b in MUTUALLY_EXCLUSIVE:
            hits_a = top_genes & set_a
            hits_b = top_genes & set_b
            # Check if both sides have high expression signals
            threshold = 0.50 if set_a == GLIA_CONFLICT_MARKERS else 0.40
            strong_a = [g for g in hits_a if gene_ratio.get(g, 0) > threshold]
            strong_b = [g for g in hits_b if gene_ratio.get(g, 0) > threshold]
            if strong_a and strong_b:
                # Determine which lineages are conflicting
                lineage_a = _infer_lineage(set_a)
                lineage_b = _infer_lineage(set_b)
                conflicts.append(f"{lineage_a}+{lineage_b}")

        if conflicts:
            formatted = []
            for conflict in conflicts:
                parts = [_format_lineage(p) for p in conflict.split("+")]
                formatted.append("+".join(parts))
            label = "/".join(formatted)
            if any("T cell+B cell" in conflict or "B cell+T cell" in conflict for conflict in formatted):
                results[cid] = f"Doublet({label})"
            else:
                results[cid] = f"Mixed({label})"
        else:
            results[cid] = "clean"
    return results


def _infer_lineage(gene_set: set[str]) -> str:
    """Map a gene set to a broad lineage name."""
    checks = [
        ("neuron", {"SNAP25", "SYT1", "RBFOX3", "SLC17A7", "GAD1", "GAD2"}),
        ("glia", {"MBP", "PLP1", "MOG", "CNP", "GFAP", "AQP4"}),
        ("T-cell", {"CD3E", "CD3D", "CD3G"}),
        ("B-cell", {"CD19", "MS4A1", "CD79A"}),
        ("endothelial", {"PECAM1", "CDH5", "VWF"}),
        ("fibroblast", {"COL1A1", "COL1A2", "DCN", "LUM"}),
        ("smooth_muscle", {"ACTA2", "MYH11", "TAGLN"}),
    ]
    for name, markers in checks:
        if gene_set & markers:
            return name
    return "unknown"


def _format_lineage(name: str) -> str:
    return CLASSIC_NAME_ALIASES.get(name, name)


def flag_ambient_rna(stats: dict, tissue: str | None = None) -> set[str]:
    """
    Identify genes that appear in every cluster at moderate levels —
    likely ambient RNA contamination, not true markers.
    """
    if not stats:
        return set()
    all_clusters = list(stats.keys())
    n_clusters = len(all_clusters)
    candidate_genes = AMBIENT_PRONE

    ambient = set()
    for gene in candidate_genes:
        cluster_count = 0
        for cid in all_clusters:
            info = stats[cid]["gene_info"].get(gene, {})
            if info and info.get("mean_expression", 0) > 0.01:
                cluster_count += 1
        # If present in >80% of clusters at moderate level → ambient
        if cluster_count > n_clusters * 0.8:
            ambient.add(gene)
    return ambient


# ═══════════════════════════════════════════════════════════════
# Dimension 4: Hierarchical Annotation
# ═══════════════════════════════════════════════════════════════

def classify_broad_class(
    top_genes: list[str],
    gene_info: dict,
    ambient_genes: set[str],
    stats_by_cluster: dict,
    cluster_id: int,
) -> tuple[str, dict[str, float]]:
    """
    Assign broad class using weighted marker hits.
    Each gene's weight = cell_ratio * specificity_factor (capped at 1.0).

    This solves the ambient RNA problem: MBP with cell_ratio=0.27 contributes
    only 0.27 weight, while DRD2 with cell_ratio=0.94 contributes 0.94.
    """
    top15 = top_genes[:15]
    # Build weighted gene scores
    gs = stats_by_cluster.get(cluster_id, {})
    gi = gs.get("gene_info", gene_info)

    def _gene_weight(g: str) -> float:
        """Weight by cell_ratio × specificity (capped at 1.0 per gene)."""
        if g not in gi:
            return 0.0
        ratio = gi[g]["cell_ratio"]
        spec = min(gi[g]["specificity"], 5.0) / 5.0  # normalize to 0-1
        return ratio * spec

    def _class_weighted_score(marker_sets: dict[str, set[str]]) -> float:
        """Sum weighted hits across all subtypes in a class."""
        all_markers = set()
        for mk in marker_sets.values():
            all_markers.update(mk)
        # Remove ambient-prone genes from consideration
        clean_markers = all_markers - ambient_genes
        total = sum(_gene_weight(g) for g in top15 if g in clean_markers)
        return round(total, 3)

    def _subtype_weighted_hits(marker_sets: dict[str, set[str]]) -> dict[str, float]:
        """Per-subtype weighted scores."""
        hits = {}
        for name, markers in marker_sets.items():
            clean_markers = markers - ambient_genes
            hits[name] = round(sum(_gene_weight(g) for g in top15 if g in clean_markers), 3)
        return hits

    # Weighted scores per broad class
    glia_score = _class_weighted_score(BROAD_CLASSES["Glia"])
    neuron_score = _class_weighted_score(BROAD_CLASSES["Neuron"])
    immune_score = _class_weighted_score(BROAD_CLASSES["Immune"])
    vascular_score = _class_weighted_score(BROAD_CLASSES["Vascular"])
    progenitor_score = _class_weighted_score(BROAD_CLASSES["Progenitor"])

    # Glia lineage markers (definitive transcription factors, less ambient-prone)
    glia_lineage = {"SOX10", "OLIG1", "OLIG2", "GFAP", "AQP4", "TMEM119", "P2RY12"}
    glia_lineage_score = sum(_gene_weight(g) for g in top15 if g in glia_lineage)

    # Neuron lineage markers
    neuron_lineage = {"SNAP25", "SYT1", "RBFOX3", "GAD1", "GAD2", "SLC17A7", "SLC17A6"}
    neuron_lineage_score = sum(_gene_weight(g) for g in top15 if g in neuron_lineage)

    scores = {
        "Glia": glia_score,
        "Neuron": neuron_score,
        "Immune": immune_score,
        "Vascular": vascular_score,
        "Progenitor": progenitor_score,
    }

    # Glia prioritization: only if definitive lineage TFs dominate
    # AND neuron lineage markers are absent/weak
    if glia_lineage_score > 1.0 and neuron_lineage_score < 0.5:
        return "Glia", scores
    if glia_lineage_score > 0.5 and neuron_lineage_score == 0:
        return "Glia", scores

    best_class = max(scores, key=scores.get)
    if scores[best_class] == 0:
        return "Unclassified", scores
    return best_class, scores


def assign_subtype(
    broad_class: str,
    top_genes: list[str],
    gene_info: dict,
    ambient_genes: set[str],
) -> tuple[str, float]:
    """
    Assign subtype using weighted gene hits (cell_ratio × specificity).
    Returns (subtype_name, confidence).
    """
    top15 = top_genes[:15]

    def _gene_weight(g: str) -> float:
        if g not in gene_info:
            return 0.0
        ratio = gene_info[g]["cell_ratio"]
        spec = min(gene_info[g]["specificity"], 5.0) / 5.0
        return ratio * spec

    def _score_subtype(markers: set[str]) -> float:
        clean = markers - ambient_genes
        return round(sum(_gene_weight(g) for g in top15 if g in clean), 3)

    if broad_class == "Glia":
        subtypes = [
            ("Oligodendrocyte", {"MBP", "PLP1", "MOG", "SOX10", "CNP", "ASPA", "CLDN11", "MAG"}),
            ("OPC", {"PDGFRA", "OLIG1", "OLIG2", "CSPG4", "GPR17", "SOX10"}),
            ("Astrocyte", {"AQP4", "GFAP", "SLC1A3", "ALDH1L1", "GJA1", "GLUL"}),
            ("Microglia", {"TMEM119", "P2RY12", "CX3CR1", "CSF1R", "ITGAM", "CD68"}),
        ]
    elif broad_class == "Neuron":
        subtypes = [
            ("D1-MSN", {"DRD1", "TAC1", "PDYN", "PPP1R1B", "ISL1", "GPR149",
                        "OTOF", "CPNE4", "SPHKAP", "RELN"}),
            ("D2-MSN", {"DRD2", "PENK", "GPR6", "ADORA2A", "GPR52"}),
            ("PV-interneuron", {"PVALB", "KCNS3", "KCNC1", "KCNC2"}),
            ("SST-interneuron", {"SST", "NPY", "CALB2", "LHX6", "CXCL14"}),
            ("VIP-interneuron", {"VIP", "CCK", "HTR3A", "CALB2"}),
            ("Cholinergic-interneuron", {"CHAT", "SLC5A7", "ACHE"}),
            ("Excitatory-neuron", {"SLC17A7", "SLC17A6", "CAMK2A", "SATB2"}),
            ("Inhibitory-neuron", {"GAD1", "GAD2", "SLC32A1"}),
            ("Dopaminergic-neuron", {"TH", "SLC6A3", "NR4A2"}),
            ("Neuroblast", {"DCX", "SOX4", "SOX11", "STMN1", "ELAVL4", "TUBB3"}),
        ]
    elif broad_class == "Immune":
        subtypes = [
            ("T-cell", {"CD3E", "CD3D", "CD3G", "CD2"}),
            ("CD4+ T-cell", {"CD4", "IL7R"}),
            ("CD8+ T-cell", {"CD8A", "CD8B", "NKG7"}),
            ("B-cell", {"CD19", "MS4A1", "CD79A", "CD79B", "PAX5"}),
            ("Plasma-cell", {"SDC1", "MZB1", "JCHAIN"}),
            ("NK-cell", {"NKG7", "GNLY", "KLRD1", "KLRF1"}),
            ("Macrophage", {"CD68", "CD163", "CSF1R"}),
            ("Monocyte", {"CD14", "S100A8", "S100A9", "LYZ"}),
            ("Dendritic-cell", {"CLEC10A", "FCER1A", "CLEC9A"}),
            ("Neutrophil", {"ELANE", "MPO", "CSF3R"}),
            ("Microglia", {"TMEM119", "P2RY12", "CX3CR1"}),
        ]
    elif broad_class == "Vascular":
        subtypes = [
            ("Endothelial", {"PECAM1", "CDH5", "VWF", "ESAM", "CLDN5"}),
            ("Pericyte", {"PDGFRB", "RGS5", "CSPG4"}),
            ("VSMC", {"ACTA2", "MYH11", "TAGLN", "CNN1"}),
        ]
    elif broad_class == "Progenitor":
        subtypes = [
            ("Proliferating", {"MKI67", "TOP2A", "PCNA", "CCNB1", "BIRC5"}),
            ("Neuroblast", {"DCX", "SOX4", "SOX11"}),
            ("NPC", {"SOX2", "NES", "PAX6", "HES5"}),
        ]
    else:
        return "Unclassified", 0.0

    best = ("Unclassified", 0.0)
    for name, markers in subtypes:
        s = _score_subtype(markers)
        if s > best[1]:
            best = (name, s)

    confidence = round(best[1], 2)
    if confidence < 0.15:
        return "Unclassified", confidence
    return best[0], confidence


def assign_reference_type(
    top_genes: list[str],
    gene_info: dict,
    ambient_genes: set[str],
    reference_markers: dict[str, set[str]],
) -> tuple[str, float, set[str]]:
    """Fallback annotation against bundled/reference marker panels."""
    top15 = top_genes[:15]

    def _gene_weight(g: str) -> float:
        if g not in gene_info:
            return 0.0
        ratio = gene_info[g]["cell_ratio"]
        spec = min(gene_info[g]["specificity"], 5.0) / 5.0
        return ratio * spec

    best_type = "Unclassified"
    best_score = 0.0
    best_markers = set()
    for cell_type, markers in reference_markers.items():
        clean = markers - ambient_genes
        score = sum(_gene_weight(g) for g in top15 if g in clean)
        if score > best_score:
            best_type = cell_type
            best_score = score
            best_markers = clean

    confidence = round(best_score, 2)
    if confidence < 0.30:
        return "Unclassified", confidence, set()
    return best_type, confidence, best_markers


# ═══════════════════════════════════════════════════════════════
# Final Score Computation
# ═══════════════════════════════════════════════════════════════

def compute_final_score(
    cluster_id: int,
    stats: dict,
    broad_class: str,
    subtype: str,
    subtype_conf: float,
    doublet_status: str,
    ambient_genes: set[str],
    canonical_markers: set[str],
) -> float:
    """
    Score = statistical + marker match + lineage purity + marker strength.
    Database concordance is not awarded by this heuristic script.
    """
    genes_ranked = stats[cluster_id]["genes_ranked"]
    gene_info = stats[cluster_id]["gene_info"]
    top10 = genes_ranked[:10]
    top15 = set(genes_ranked[:15])

    # Statistical support (0.30): top marker specificity should not be diluted
    # by long-tail background genes later in the marker table.
    top_stat_genes = top10[:5]
    passing = sum(1 for g in top_stat_genes if gene_info[g]["passes_ratio_filter"])
    stat_score = (passing / max(len(top_stat_genes), 1)) * 0.30

    # Marker match (0.30): canonical markers found in top 10
    if canonical_markers:
        hits = len(set(top10) & canonical_markers)
        if hits >= 3:
            marker_score = 0.30
        elif hits == 2:
            marker_score = 0.20
        elif hits == 1:
            marker_score = 0.10
        else:
            marker_score = 0.0
    else:
        marker_score = 0.05  # no markers available, minimal credit

    # Lineage purity (0.25): check for contradictory lineage markers, minus ambient
    clean_top15 = top15 - ambient_genes
    contradictory = 0
    if broad_class == "Glia":
        conflict_set = (BROAD_CLASSES["Neuron"]["generic"] |
                        BROAD_CLASSES["Neuron"]["excitatory"] |
                        BROAD_CLASSES["Neuron"]["inhibitory"] |
                        BROAD_CLASSES["Immune"]["t_cell"] |
                        BROAD_CLASSES["Immune"]["b_cell"])
        contradictory = len(clean_top15 & conflict_set)
    elif broad_class == "Neuron":
        conflict_set = (BROAD_CLASSES["Glia"]["oligodendrocyte"] |
                        BROAD_CLASSES["Glia"]["astrocyte"] |
                        BROAD_CLASSES["Glia"]["microglia"] |
                        BROAD_CLASSES["Vascular"]["endothelial"])
        contradictory = len(clean_top15 & conflict_set)
    elif broad_class == "Immune":
        conflict_set = (BROAD_CLASSES["Neuron"]["generic"] |
                        BROAD_CLASSES["Glia"]["oligodendrocyte"])
        contradictory = len(clean_top15 & conflict_set)

    if contradictory == 0:
        purity_score = 0.25
    elif contradictory <= 2:
        purity_score = 0.12
    else:
        purity_score = 0.0

    # Doublet penalty
    if doublet_status != "clean":
        purity_score = min(purity_score, 0.05)

    # Marker strength (0.15): subtype-specific weighted signal.
    marker_strength_score = min(max(subtype_conf, 0.0), 2.0) / 2.0 * 0.15

    # Full database/literature concordance requires external validation.
    db_score = 0.0

    total = round(stat_score + marker_score + purity_score + marker_strength_score + db_score, 2)
    return min(total, 1.0)


def default_output_path(input_path: str) -> str:
    """Create a non-destructive default output path next to the input file."""
    path = Path(input_path)
    return str(path.with_name("anno_candidate_annotation.txt"))


def confidence_label(score: float) -> str:
    if score >= 0.80:
        return "high"
    if score >= 0.50:
        return "moderate"
    if score >= 0.30:
        return "low"
    return "unreliable"


def _base_annotation_name(label: str) -> str:
    base = label.split("(", 1)[0]
    base = re.sub(r"-\d+$", "", base)
    return base


def _marker_panel_for_label(base: str) -> set[str]:
    if base == "Endothelial":
        base = "Endothelial cell"
    if base in BRAIN_TYPE_MARKERS:
        return BRAIN_TYPE_MARKERS[base]
    if "interneuron" in base.lower():
        return BRAIN_TYPE_MARKERS["Interneuron"]
    return set()


def marker_basis_for(cluster_id: int, stats: dict, max_markers: int = 5) -> str:
    markers = stats[cluster_id]["genes_ranked"][:max_markers]
    return "+".join(markers)


def _is_brain_tissue(tissue: str | None) -> bool:
    tissue_norm = (tissue or "").lower().replace("-", "_").replace(" ", "_")
    return tissue_norm in {
        "brain", "brain_multi_region", "multi_region_brain",
        "cortex_subcortex", "striatum", "cortex",
    }


def _gene_weight_from_stats(
    cluster_id: int,
    stats: dict,
    gene: str,
    ambient_genes: set[str] | None = None,
) -> float:
    if ambient_genes and gene in ambient_genes:
        return 0.0
    info = stats[cluster_id]["gene_info"].get(gene, {})
    if not info:
        return 0.0
    ratio = info.get("cell_ratio", 0.0)
    specificity = min(info.get("specificity", 0.0), 5.0) / 5.0
    return ratio * specificity


def _panel_signal(
    cluster_id: int,
    stats: dict,
    markers: set[str],
    ambient_genes: set[str] | None = None,
    top_n: int = 15,
) -> float:
    top = set(stats[cluster_id]["genes_ranked"][:top_n])
    return round(sum(
        _gene_weight_from_stats(cluster_id, stats, gene, ambient_genes)
        for gene in top & markers
    ), 3)


def _strong_marker_hits(
    cluster_id: int,
    stats: dict,
    markers: set[str],
    min_ratio: float = 0.30,
    max_hits: int = 4,
) -> list[str]:
    hits = []
    for gene in stats[cluster_id]["genes_ranked"][:15]:
        if gene not in markers:
            continue
        info = stats[cluster_id]["gene_info"].get(gene, {})
        if info.get("cell_ratio", 0.0) >= min_ratio:
            hits.append(gene)
    return hits[:max_hits]


def _best_panel_label(cluster_id: int, stats: dict, labels: list[str], ambient_genes: set[str]) -> tuple[str, float]:
    scored = [
        (label, _panel_signal(cluster_id, stats, BRAIN_TYPE_MARKERS.get(label, set()), ambient_genes))
        for label in labels
    ]
    return max(scored, key=lambda item: item[1])


def apply_brain_conservative_rules(
    subtype: str,
    broad_class: str,
    cluster_id: int,
    stats: dict,
    ambient_genes: set[str],
    tissue: str | None,
) -> tuple[str, str, float]:
    """Prefer lineage-clean glia; flag balanced glia-neuron signal as mixed."""
    if not _is_brain_tissue(tissue):
        return subtype, broad_class, 0.0

    glia_labels = ["Oligodendrocyte", "OPC", "Astrocyte", "Microglia"]
    neuron_labels = ["D1-MSN", "D2-MSN", "Interneuron", "Excitatory-neuron", "Neuroblast"]
    best_glia, glia_signal = _best_panel_label(cluster_id, stats, glia_labels, ambient_genes)
    best_neuron, neuron_signal = _best_panel_label(cluster_id, stats, neuron_labels, ambient_genes)

    glia_hits = _strong_marker_hits(cluster_id, stats, BRAIN_TYPE_MARKERS[best_glia])
    neuron_hits = _strong_marker_hits(cluster_id, stats, BRAIN_TYPE_MARKERS[best_neuron])

    if glia_signal >= 0.90 and neuron_signal >= 0.75 and glia_hits and neuron_hits:
        ratio = min(glia_signal, neuron_signal) / max(glia_signal, neuron_signal)
        if ratio >= 0.55:
            return "Mixed(Glia+Neuron)", "Mixed", min(glia_signal, neuron_signal)

    if glia_signal >= 0.55 and glia_signal >= neuron_signal * 1.15 and glia_hits:
        return best_glia, "Glia", glia_signal

    if broad_class == "Neuron" and subtype in BRAIN_TYPE_MARKERS:
        subtype_signal = _panel_signal(cluster_id, stats, BRAIN_TYPE_MARKERS[subtype], ambient_genes)
        if glia_signal >= 0.55 and glia_signal > subtype_signal and glia_hits:
            return best_glia, "Glia", glia_signal

    if subtype != "Unclassified" and subtype in BRAIN_TYPE_MARKERS:
        return subtype, broad_class, _panel_signal(
            cluster_id, stats, BRAIN_TYPE_MARKERS[subtype], ambient_genes,
        )
    return subtype, broad_class, 0.0


def evidence_support_for(label: str, cluster_id: int, stats: dict) -> str:
    base = _base_annotation_name(label)
    if base == "Endothelial":
        base = "Endothelial cell"
    markers = _marker_panel_for_label(base)
    hits = _strong_marker_hits(cluster_id, stats, markers, min_ratio=0.25)
    if "Doublet" in base:
        return "Candidate: lineage-conflict rule; verify with marker/literature review"
    if "Mixed" in base:
        return "Candidate: conflicting lineage marker sets; QC/doublet review and database/literature validation required"
    if hits:
        hit_text = "/".join(hits)
        if base in {"D1-MSN", "D2-MSN", "Interneuron", "Excitatory-neuron", "Neuroblast"} or "interneuron" in base.lower():
            return f"Candidate: {hit_text} brain-neuron marker panel; validate with Allen, CellMarker 2.0, PanglaoDB, PubMed"
        if base in {"Oligodendrocyte", "OPC", "Astrocyte", "Microglia"}:
            return f"Candidate: {hit_text} glial marker panel; validate with CellMarker 2.0, PanglaoDB, Allen/PubMed"
        if base in {"Endothelial cell", "Vascular smooth muscle cell", "Pericyte"}:
            return f"Candidate: {hit_text} vascular marker panel; validate with CellMarker 2.0/PanglaoDB"
        if base in {"T-cell", "B-cell", "NK-cell"}:
            return f"Candidate: {hit_text} immune marker panel; validate with CellMarker 2.0/PanglaoDB"
    return "Candidate: marker-panel heuristic; PubMed/database review required"


def annotation_note_for(label: str, score: float, cluster_id: int, stats: dict) -> str:
    base = _base_annotation_name(label)
    if base == "Endothelial":
        base = "Endothelial cell"
    markers = set(stats[cluster_id]["genes_ranked"][:8])
    glia_signal = max(
        _panel_signal(cluster_id, stats, BRAIN_TYPE_MARKERS[label], set())
        for label in ["Oligodendrocyte", "OPC", "Astrocyte", "Microglia"]
    )
    neuron_signal = max(
        _panel_signal(cluster_id, stats, BRAIN_TYPE_MARKERS[label], set())
        for label in ["D1-MSN", "D2-MSN", "Interneuron", "Excitatory-neuron", "Neuroblast"]
    )
    if score < 0.30:
        return "Very low score; treat as provisional and inspect marker specificity"
    if score < 0.50:
        return "Low score; use as candidate annotation pending manual/literature review"
    if "Doublet" in base or "Mixed" in base:
        return "Strong markers from different lineages co-occur; review doublet, ambient RNA, and UMAP location"
    if base in {"D1-MSN", "D2-MSN", "Interneuron", "Excitatory-neuron"} and glia_signal >= neuron_signal * 0.6:
        return "Neuron-like call with substantial glial signal; keep provisional until QC and literature/database review"
    if base in {"Oligodendrocyte", "OPC", "Astrocyte", "Microglia"} and glia_signal > neuron_signal:
        return "Glial lineage markers dominate; weaker neuronal markers should be treated as possible background"
    if base == "D1-MSN":
        return "MSN-like cluster; D1 assignment supported by D1/direct-pathway markers when present"
    if base == "D2-MSN":
        return "MSN-like cluster; D2 assignment supported by PENK/DRD2/GPR6/ADORA2A when present"
    if base == "Oligodendrocyte":
        return "Oligodendrocyte-lineage markers dominate; check ambient myelin markers across clusters"
    if base == "OPC":
        return "OPC call supported by PDGFRA/OLIG1/OLIG2/GPR17-like markers"
    if base == "Astrocyte":
        return "Astrocyte call supported by AQP4/SLC1A3/GFAP/ALDH1L1-like markers"
    if base == "Microglia":
        return "Microglia call supported by TMEM119/P2RY12/CX3CR1-like markers"
    if base == "Neuroblast":
        return "Immature/progenitor-like neuronal markers; inspect DCX/SOX4/SOX11/cell-cycle signal"
    if base == "Endothelial cell":
        return "Vascular endothelial markers dominate"
    if markers & {"SST", "PVALB", "VIP", "LHX6", "NPY"}:
        return "Interneuron-like marker panel; subtype should be checked against regional context"
    return "Annotation based on marker specificity, tissue prior, and local marker panel"


def _canonicalize_cell_name(cell_type: str) -> str:
    return CLASSIC_NAME_ALIASES.get(cell_type, cell_type)


def _append_qualifier(base_name: str, qualifier: str) -> str:
    if not qualifier:
        return base_name
    if base_name.endswith(")") and "(" in base_name:
        return f"{base_name[:-1]}; {qualifier})"
    return f"{base_name}({qualifier})"


def _detect_state_qualifier(cluster_id: int, stats: dict, base_name: str) -> str:
    top15 = set(stats[cluster_id]["genes_ranked"][:15])
    base_norm = base_name.lower()
    for state, markers in STATE_SIGNATURES:
        if state in base_norm:
            continue
        if len(top15 & markers) >= 2:
            return state
    return ""


def _candidate_markers_for_duplicate(
    cluster_id: int,
    sibling_ids: list[int],
    stats: dict,
    base_type: str = "",
) -> list[str]:
    preferred_panel = _marker_panel_for_label(_base_annotation_name(base_type))
    sibling_presence = defaultdict(int)
    sibling_top = {}
    for sid in sibling_ids:
        genes = stats[sid]["genes_ranked"][:25]
        sibling_top[sid] = genes
        for gene in set(genes):
            sibling_presence[gene] += 1

    def _rank_candidates(allowed: set[str]) -> list[tuple[float, str]]:
        ranked = []
        gene_info = stats[cluster_id]["gene_info"]
        for gene in sibling_top[cluster_id]:
            if allowed and gene not in allowed:
                continue
            if gene in LOW_VALUE_QUALIFIER_MARKERS:
                continue
            info = gene_info.get(gene, {})
            if not info:
                continue
            if info.get("cell_ratio", 0.0) < 0.50:
                continue
            specificity = info.get("specificity", 0.0)
            unique_bonus = 2.0 if sibling_presence[gene] == 1 else 1.0
            if allowed or sibling_presence[gene] == 1 or specificity >= 2.0:
                score = info.get("rank_score", 0.0) * unique_bonus
                ranked.append((score, gene))
        return ranked

    ranked = _rank_candidates(preferred_panel)
    if not ranked:
        ranked = _rank_candidates(set())

    ranked.sort(reverse=True)
    seen = set()
    markers = []
    for _, gene in ranked:
        if gene in seen:
            continue
        markers.append(gene)
        seen.add(gene)
    return markers


def _fallback_markers_for_name(cluster_id: int, stats: dict) -> list[str]:
    markers = [
        gene for gene in stats[cluster_id]["genes_ranked"][:10]
        if gene not in LOW_VALUE_QUALIFIER_MARKERS
    ]
    if markers:
        return markers
    return stats[cluster_id]["genes_ranked"][:3]


def infer_candidate_label_prior(
    cluster_id: int,
    cluster_rows: list[dict],
    stats: dict,
) -> tuple[str | None, float, float]:
    """Use the input cell_type column as a weak candidate-label prior."""
    ranked_genes = stats[cluster_id]["genes_ranked"]
    rank_index = {gene: idx for idx, gene in enumerate(ranked_genes)}
    scores = defaultdict(float)

    for row in cluster_rows:
        raw_label = row.get("cell_type", "")
        if not raw_label or raw_label in {"unknown", "FindAllMarkers"}:
            continue
        label = _canonicalize_cell_name(raw_label)
        gene = row["gene"]
        info = stats[cluster_id]["gene_info"].get(gene, {})
        if not info:
            continue
        rank = rank_index.get(gene, 999)
        rank_weight = 1.0 + (1.0 / ((rank + 1) ** 0.5))
        specificity = min(info.get("specificity", 0.0), 5.0) / 5.0
        expression = min(info.get("mean_expression", 0.0), 3.0) / 3.0
        ratio = info.get("cell_ratio", 0.0)
        weight = ratio * max(specificity, 0.15) * (1.0 + expression) * rank_weight
        if ratio > 0.60:
            weight *= 1.25
        scores[label] += weight

    if not scores:
        return None, 0.0, 0.0
    ordered = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    top_label, top_score = ordered[0]
    second_score = ordered[1][1] if len(ordered) > 1 else 0.0
    return top_label, round(top_score, 3), round(second_score, 3)


def apply_candidate_label_prior(
    subtype: str,
    broad_class: str,
    cluster_id: int,
    cluster_rows: list[dict],
    stats: dict,
) -> tuple[str, str]:
    candidate, cand_score, second_score = infer_candidate_label_prior(
        cluster_id, cluster_rows, stats,
    )
    if not candidate:
        return subtype, broad_class

    candidate_broad = CANDIDATE_LABEL_BROAD_CLASS.get(candidate)
    if not candidate_broad:
        return subtype, broad_class

    top_gene = stats[cluster_id]["genes_ranked"][0] if stats[cluster_id]["genes_ranked"] else ""
    top_gene_label = None
    for row in cluster_rows:
        if row["gene"] == top_gene:
            top_gene_label = _canonicalize_cell_name(row.get("cell_type", ""))
            break

    subtype_norm = _canonicalize_cell_name(subtype)
    same_broad = candidate_broad == broad_class or broad_class in {"Unclassified", "Reference"}
    strong_prior = cand_score >= 1.0 and cand_score >= max(second_score * 1.10, 0.6)
    top_gene_ratio = stats[cluster_id]["gene_info"].get(
        top_gene, {},
    ).get("cell_ratio", 0.0)
    top_gene_support = top_gene_label == candidate and top_gene_ratio > 0.60
    top_gene_broad = CANDIDATE_LABEL_BROAD_CLASS.get(top_gene_label or "")
    if top_gene_label and top_gene_broad == broad_class and top_gene_ratio > 0.75:
        return top_gene_label, top_gene_broad

    if subtype_norm == "Unclassified" and cand_score >= 0.15:
        return candidate, candidate_broad
    if same_broad and (strong_prior or top_gene_support):
        return candidate, candidate_broad
    return subtype, broad_class


def make_unique_annotation_names(
    cluster_types: dict[int, str],
    stats: dict,
) -> dict[int, str]:
    """
    Return unique, concise cluster labels.
    Keep the classic base name when it is already unique. Add gene suffixes
    only to distinguish duplicated base names. If marker suffixes still collide,
    add -1/-2 after the base name.
    """
    type_counts = defaultdict(list)
    for cid, ctype in cluster_types.items():
        type_counts[_canonicalize_cell_name(ctype)].append(cid)

    final_names = {}
    for ctype, cids in type_counts.items():
        cids_sorted = sorted(cids)
        if len(cids_sorted) == 1:
            cid = cids_sorted[0]
            final_names[cid] = ctype
            continue

        initial_candidates = {}
        candidate_counts = defaultdict(int)
        for cid in cids_sorted:
            markers = _candidate_markers_for_duplicate(cid, cids_sorted, stats, ctype)
            if not markers and not _marker_panel_for_label(_base_annotation_name(ctype)):
                markers = _fallback_markers_for_name(cid, stats)
            suffix = "+".join(sorted(markers[:2])) if markers else ""
            candidate = _append_qualifier(ctype, suffix) if suffix else ctype
            initial_candidates[cid] = (candidate, suffix)
            candidate_counts[candidate] += 1

        needs_ordinals = (
            any(count > 1 for count in candidate_counts.values()) or
            any(not suffix for _, suffix in initial_candidates.values())
        )
        for ordinal, cid in enumerate(cids_sorted, start=1):
            candidate, suffix = initial_candidates[cid]
            if not needs_ordinals:
                final_names[cid] = candidate
                continue
            if suffix:
                final_names[cid] = f"{ctype}-{ordinal}({suffix})"
            else:
                final_names[cid] = f"{ctype}-{ordinal}"

    return final_names


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Four-dimension scRNA-seq annotation scoring")
    parser.add_argument("input", help="Input TSV (FindAllMarkers format)")
    parser.add_argument("--tissue", default=None,
                        help="Tissue context (auto-detect if omitted)")
    parser.add_argument("--species", default="human", choices=["human", "mouse"])
    parser.add_argument("--output", default=None, help="Output file path")
    parser.add_argument("--markers", help="Additional marker TSV")
    parser.add_argument("--markers-json", help="Additional marker JSON")
    args = parser.parse_args()

    outpath = args.output or default_output_path(args.input)
    if Path(outpath).resolve() == Path(args.input).resolve():
        raise ValueError("Output path resolves to the input file; refusing to overwrite input")

    rows, input_header = load_input(args.input, species=args.species)

    # Load additional markers if provided
    extra_markers = {}
    if args.markers_json:
        with open(args.markers_json) as f:
            extra_markers = {k: set(v) for k, v in json.load(f).items()}
    elif args.markers:
        extra_markers = load_markers(args.markers)

    # Load bundled reference
    script_dir = Path(__file__).parent.parent
    markers_path = script_dir / "references" / "markers.md"
    bundled_markers = {}
    if markers_path.exists():
        bundled_markers = parse_markers_md(str(markers_path))

    # Merge
    all_markers = {**bundled_markers, **extra_markers}

    # Dimension 1: statistics
    stats = compute_cluster_stats(rows)

    # Dimension 2: tissue inference
    tissue = args.tissue or infer_tissue(rows)
    print(f"# Inferred tissue: {tissue}")

    # Dimension 3: ambient RNA and doublets
    cluster_data_raw = defaultdict(list)
    for row in rows:
        cluster_data_raw[row["cluster"]].append(row)

    ambient_genes = flag_ambient_rna(stats, tissue)
    if ambient_genes:
        print(f"# Flagged ambient RNA genes: {', '.join(sorted(ambient_genes))}")

    cluster_genes_for_doublets = {
        cid: data["genes_ranked"] for cid, data in stats.items()
    }
    doublet_results = detect_doublets(
        cluster_genes_for_doublets, cluster_data_raw, ambient_genes,
    )

    # Group by cluster
    cluster_genes = defaultdict(list)
    cluster_prelim = {}
    for row in rows:
        cid = row["cluster"]
        cluster_genes[cid].append(row["gene"])
        if cid not in cluster_prelim:
            cluster_prelim[cid] = row["cell_type"]

    # Dimension 4: hierarchical annotation
    cluster_types = {}
    cluster_scores = {}
    cluster_classes = {}

    for cid in sorted(cluster_genes.keys()):
        genes_sorted = stats[cid]["genes_ranked"]
        gene_info = stats[cid]["gene_info"]

        # Step 4a: Broad class (weighted by statistical quality, ambient-excluded)
        broad_class, class_scores = classify_broad_class(
            genes_sorted, gene_info, ambient_genes, stats, cid,
        )

        # Step 4b: Subtype (weighted)
        subtype, subtype_conf = assign_subtype(
            broad_class, genes_sorted, gene_info, ambient_genes,
        )

        reference_markers = set()
        if subtype == "Unclassified":
            ref_type, ref_conf, reference_markers = assign_reference_type(
                genes_sorted, gene_info, ambient_genes, all_markers,
            )
            if ref_type != "Unclassified":
                subtype = ref_type
                subtype_conf = ref_conf
                if broad_class == "Unclassified":
                    broad_class = "Reference"

        subtype, broad_class = apply_candidate_label_prior(
            subtype, broad_class, cid, cluster_data_raw.get(cid, []), stats,
        )

        subtype, broad_class, rule_conf = apply_brain_conservative_rules(
            subtype, broad_class, cid, stats, ambient_genes, tissue,
        )
        subtype_conf = max(subtype_conf, rule_conf)

        if subtype.startswith("Mixed("):
            cluster_types[cid] = subtype
            cluster_scores[cid] = 0.32
            cluster_classes[cid] = "Mixed"
            continue

        # Handle doublets
        if doublet_results.get(cid, "clean") != "clean":
            cluster_types[cid] = doublet_results[cid]
            cluster_scores[cid] = 0.3
            cluster_classes[cid] = "Mixed"
            continue

        # Collect canonical markers for this type
        canonical = set()
        for ct, mk in all_markers.items():
            ct_norm = ct.lower().replace(" ", "").replace("-", "")
            st_norm = subtype.lower().replace(" ", "").replace("-", "")
            if st_norm in ct_norm or ct_norm in st_norm:
                canonical.update(mk)
        canonical.update(reference_markers)
        # Also add subtype markers from built-in definitions
        canonical.update(_get_subtype_markers(subtype))

        # Final score
        score = compute_final_score(
            cid, stats, broad_class, subtype, subtype_conf,
            doublet_results.get(cid, "clean"), ambient_genes, canonical,
        )

        cluster_classes[cid] = broad_class
        cluster_types[cid] = subtype
        cluster_scores[cid] = score

    # Ensure each cluster receives a unique, interpretable label.
    final_names = make_unique_annotation_names(cluster_types, stats)

    # Output
    annotation_cols = [
        "manual_type", "score", "confidence",
        "marker_basis", "evidence_support", "annotation_note",
    ]
    stale_cols = {
        # Historical typo is removed from old inputs; never emitted.
        "manual_type", "mannual_type", "score", "confidence",
        "marker_basis", "evidence_support", "annotation_note",
        "broad_class", "evidence_status", "naming_status",
    }
    output_header = [col for col in input_header if col not in stale_cols] + annotation_cols
    with open(outpath, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(output_header)
        for row in rows:
            cid = row["cluster"]
            raw = row.get("_raw", {})
            score = cluster_scores.get(cid, 0.0)
            label = final_names.get(cid, "Unclassified")
            additions = {
                "manual_type": label,
                "score": str(score),
                "confidence": confidence_label(score),
                "marker_basis": marker_basis_for(cid, stats),
                "evidence_support": evidence_support_for(label, cid, stats),
                "annotation_note": annotation_note_for(label, score, cid, stats),
            }
            writer.writerow([
                additions[col] if col in additions else raw.get(col, "")
                for col in output_header
            ])

    # Print summary
    print(f"\n## Annotation Summary ({len(cluster_genes)} clusters, tissue: {tissue}, species: {args.species})")
    print(f"{'Cl':<4} {'Original':<22} {'Broad':<14} {'manual_type':<32} {'Score':<6} {'Conf':<10} {'Top Markers'}")
    print("-" * 108)
    for cid in sorted(cluster_genes.keys()):
        top3 = ", ".join(stats[cid]["genes_ranked"][:3])
        db_flag = " [DOUBLET]" if doublet_results.get(cid, "clean") != "clean" else ""
        print(f"{cid:<4} {cluster_prelim[cid]:<22} {cluster_classes.get(cid, '?'):<14} "
              f"{final_names.get(cid, '?'):<32} {cluster_scores.get(cid, 0):<6.2f} "
              f"{confidence_label(cluster_scores.get(cid, 0)):<10} {top3}{db_flag}")

    print(f"\nOutput: {outpath}")


def _get_subtype_markers(subtype: str) -> set[str]:
    """Get built-in markers for a subtype."""
    for broad, subtypes in BROAD_CLASSES.items():
        for st, markers in subtypes.items():
            st_norm = st.lower().replace("_", "-").replace(" ", "")
            sub_norm = subtype.lower().replace("_", "-").replace(" ", "")
            if st_norm in sub_norm or sub_norm in st_norm:
                return markers
    return set()


if __name__ == "__main__":
    main()
