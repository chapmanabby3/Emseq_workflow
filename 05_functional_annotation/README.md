# Functional Annotation

This section assigns standardized gene symbols to differentially methylated loci identified in *Boreogadus saida*.

Because the reference genome annotation consists primarily of locus identifiers (e.g. `BORSAIPOL_LOCUS####`), protein homology searches were used to infer gene identities based on curated protein annotations.

## Overview

The annotation workflow consists of two steps:

1. Protein homology search against the UniProtKB/Swiss-Prot database using BLASTP.
2. Mapping BLAST results back to genome locus identifiers and assigning gene symbols to annotated DMCs.

## Workflow

### 1. Protein homology search

Predicted protein sequences from the *B. saida* genome assembly (`.faa`) are compared against the UniProtKB/Swiss-Prot protein database using BLASTP.

Swiss-Prot was selected because it contains manually reviewed protein annotations with standardized gene symbols.

For each predicted *B. saida* protein, the best-scoring homolog is retained.

### 2. Mapping proteins to genomic loci

The NCBI feature table is used to associate genome locus tags (`BORSAIPOL_LOCUS####`) with their corresponding predicted protein accession numbers (`CAL82XXXXX`).

This provides the link between:

```
Genome locus
        ↓
Predicted protein
        ↓
Best Swiss-Prot match
        ↓
Standard gene symbol
```

### Gene symbol assignment

Gene symbols are assigned only when the best BLAST hit satisfies both of the following criteria:

- Percent identity ≥ 40%
- E-value ≤ 1 × 10⁻⁵

These thresholds provide reasonable confidence that the matched protein represents a functional homolog.

When no reliable homolog is identified, the original *B. saida* locus identifier is retained.

## Inputs

- Predicted protein FASTA (`.faa`)
- NCBI feature table
- UniProtKB/Swiss-Prot protein database
- Annotated DMC table

## Outputs

The final output contains the original DMC annotation together with:

- matched protein accession
- UniProt gene symbol
- protein identity (%)
- BLAST E-value
- final gene label used for downstream analyses

These gene labels are subsequently used for ortholog mapping and functional enrichment analyses.