# Marker Gene Databases & Query Guide

Quick reference for searching and validating cell type marker genes online.

## Primary Databases

### 1. PanglaoDB
- **URL**: https://panglaodb.se/
- **Paper**: https://pubmed.ncbi.nlm.nih.gov/30951143/
- **API/Query**: `https://panglaodb.se/search?query=GENE_NAME`
- **Content**: Cell type markers from >1300 scRNA-seq studies across human and mouse
- **Best for**: Quick lookup of which cell types a gene marks
- **Usage**: Submit individual genes; returns cell type enrichments with p-values

**Web search pattern**: `site:panglaodb.se GENE` or `PanglaoDB GENE cell type`

### 2. CellMarker 2.0
- **URL**: http://bio-bigdata.hrbmu.edu.cn/CellMarker/
- **Paper**: https://academic.oup.com/nar/article/51/D1/D870/6775381
- **Content**: Manually curated cell markers from literature, covering human and mouse across 190+ tissues
- **Best for**: Authoritative, literature-backed marker validation
- **Search**: By gene, cell type, or tissue

**Web search pattern**: `CellMarker GENE cell type marker human` or search the database directly via WebFetch

### 3. Human Cell Atlas (HCA)
- **URL**: https://www.humancellatlas.org/
- **Data portal**: https://data.humancellatlas.org/
- **Best for**: Reference expression across human tissues
- **Usage**: Search for gene expression across tissues and cell types

### 4. Allen Brain Atlas / Cell Types Database
- **URL**: https://portal.brain-map.org/
- **Cell types**: https://celltypes.brain-map.org/mouse/
- **Best for**: Brain/nervous system cell type validation
- **Content**: In situ hybridization and scRNA-seq data for human and mouse brain

### 5. Cell Taxonomy (NCI)
- **URL**: https://ngdc.cncb.ac.cn/celltaxonomy/
- **Content**: Comprehensive cell type ontology with markers across species

## Supplementary Databases

### Tabula Muris / Tabula Sapiens
- **Tabula Muris**: Mouse cell atlas; markers across 20 organs
- **Tabula Sapiens**: Human cell atlas; multi-organ reference

### scMatch / CellTypist
- Reference-based automated annotation tools
- Use their published marker lists for validation

## Search Strategy per Gene

For each gene, execute searches in this order:

```
1. "GENE_NAME cell type marker" → get general consensus
2. "GENE_NAME expression cell type TISSUE" → tissue-specific check
3. "GENE_NAME marker SPECIES scRNA-seq" → validate in species
4. "canonical markers for CELL_TYPE SPECIES" → verify the full panel
```

### Web search query templates

**Single gene lookup:**
```
What cell type is GENE a marker for in SPECIES TISSUE?
```

**Cell type verification:**
```
Canonical marker genes for CELL_TYPE in SPECIES TISSUE scRNA-seq
```

**Differential expression check:**
```
Is GENE specifically expressed in CELL_TYPE or is it also expressed elsewhere?
```

**Cross-species check:**
```
GENE cell type marker human vs mouse expression
```

**Doublet check:**
```
Can GENE_A and GENE_B be co-expressed in the same cell type?
```

## When to Trust Database vs Literature

- **Prefer CellMarker** for well-studied markers (literature-curated, species-specific)
- **Prefer PanglaoDB** for novel or rarely-studied genes (broader coverage from scRNA-seq studies)
- **Prefer PubMed/literature** when databases disagree or for rare cell types
- **Cross-validate**: if 2+ databases agree on a marker assignment, confidence is higher

## Confirming Doublets / Mixed Clusters (cell-level)

Marker-table conflict is only a **candidate** doublet — it cannot distinguish a true
cell-mixing event from ambient-RNA bleed at the cluster level. Before committing to (or
dismissing) a `Doublet`/`Mixed` call, confirm with cell-resolution evidence:

- **Specificity first**: a real doublet shows *both* lineages cluster-specific (high in
  this cluster, low elsewhere). A second lineage that is high but present across most
  clusters / low-specificity is ambient RNA (e.g. striatal `PENK`, myelin `PLP1/MBP`) —
  not a doublet.
- **Per-cell doublet scores**: run **Scrublet** or **DoubletFinder** and check whether
  the cluster is enriched for high-doublet-score cells. The scoring script's `Doublet QC`
  note marks clusters that need this check.
- **Co-expression at single-cell level**: in a true doublet the two marker programs fire
  in the *same cells*; in ambient contamination the background gene is spread thinly
  across many cells. Inspect a co-expression scatter / feature plot.
- **UMAP / cluster topology**: doublet clusters frequently sit on the bridge between two
  parent clusters and may shrink or vanish after doublet removal and re-clustering.

**Web search patterns:**
```
Scrublet DoubletFinder detect doublets scRNA-seq cluster
Can CELLTYPE_A and CELLTYPE_B markers co-express same cell doublet vs ambient RNA
ambient RNA contamination SoupX correction GENE bleed scRNA-seq
```
