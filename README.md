# Permutation-Based Rescue of High-Confidence Cells in Spatial Proteomics (NSCLC)

A Python workflow that uses label-permutation testing to filter noisy cells and
rescue high-confidence cells for oncology marker analysis in a non-small cell
lung cancer (NSCLC) CODEX/PhenoCycler dataset.

## Problem
In sparse UMAP embeddings of single-cell spatial proteomics data, many cells
with real biological signal get discarded along with noisy ones. This project
tests whether a cell's local neighborhood is enriched for its assigned label
beyond what chance would produce, so valuable cells can be kept.

## Approach
1. Preprocessing: normalization / marker selection / QC steps
2. UMAP embedding and neighborhood graph construction k = X neighbors
3. Neighborhood enrichment score per cell using a matrix operation (Y = WB)
4. Label permutation test ([N] permutations) to build a null distribution
5. Scoring: median/MAD-based scores against the null
6. Empirical p-values with Benjamini-Hochberg FDR correction**
7. Rescue: keep cells significant at FDR < [0.05]

## Results
Using a small(50,000 cells) portion of the data, 
- 34,255 of 50,000 cells rescued (68.51%) compared with a simple threshold filter

## Tech stack
R, Python, NumPy, pandas, SciPy, scanpy/umap-learn, matplotlib [edit to match the notebook]

## Repository contents
- `[permutation_test.ipynb]`: used cleaned data to perform permutation analysis 
- `[umbizo seurat object]` : creates the object and UMAP matrices using raw data

## How to run
```bash
pip install -r requirements.txt
jupyter notebook [permutation_test.ipynb]
```

## Author
Jimena Roces, UCLA Computational & Systems Biology
