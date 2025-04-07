# vdsv
# CellBender + Label Transfer Pipeline for Frog Spinal Cord Neurons

This repository contains analysis code used to generate **Figure X** and **Supplementary Figure Y** of the manuscript "[TITLE]" (co-authored by Yuri). The analysis focuses on identifying spinal interneuron and motor neuron subtypes across Xenopus developmental stages, using CellBender to clean ambient RNA and Seurat for downstream processing and label transfer.

## Overview

We processed four stages: 38, 47, 54, and 57. The pipeline includes:

1. **Loading CellBender-cleaned data** (`*.h5` files) and converting to Seurat objects
2. **Filtering cells** using bimodal thresholding on `nFeature_RNA` and `nCount_RNA`
3. **Normalizing and scaling** using Seurat
4. **Label transfer**:
   - Coarse cell type annotations from stage 54 reference
   - Cardinal interneuron class transfer onto identified neurons
5. **Visualization** of UMAPs, feature plots, dot plots, and cluster proportions across replicates

Note: Some reference annotations (e.g., the stage 54 reference for coarse cell types and the cardinal neuron class reference) were generated previously and are not included in this repository.

## File Guide

- `load_cellbender_h5_to_seurat.R`: loads CellBender `.h5` output and converts it to a Seurat object
- `filter_bimodal_thresholding.R`: applies `mixtools` to determine cutoff thresholds for filtering cells
- `label_transfer_coarse_and_cardinal.R`: performs label transfer from reference datasets using Seurat anchors
- `plot_figures_main_supp.R`: generates the UMAPs and dot plots used in main and supplementary figures

## Assumptions

- This code assumes you have already run CellBender and have files named like `output_FPR_0.1_filtered.h5`
- When using your own data, assign your Seurat object to a variable named `seu` or similar for consistency across scripts

## Citation

This code accompanies the manuscript "[TITLE]" — if using this code for related work, please cite the paper when available.

Inside load_cellbender_h5_to_seurat.R
r
Copy
Edit
# Example script to load CellBender .h5 file and convert to Seurat
library(Seurat)
library(rhdf5)
library(Matrix)

# Modify path to match your CellBender output
h5file <- "/path/to/output_FPR_0.1_filtered.h5"

# Load datasets
barcodes <- h5read(h5file, "/matrix/barcodes")
data     <- h5read(h5file, "/matrix/data")
indices  <- h5read(h5file, "/matrix/indices") + 1
indptr   <- h5read(h5file, "/matrix/indptr")
shape    <- h5read(h5file, "/matrix/shape")
genes    <- h5read(h5file, "/matrix/features/name")

# Build sparse matrix
mat <- sparseMatrix(
  i = indices,
  p = indptr,
  x = data,
  dims = c(shape[1], shape[2])
)

rownames(mat) <- make.unique(genes)
colnames(mat) <- make.unique(barcodes)

# Create Seurat object
seu <- CreateSeuratObject(mat)

# You can now apply filtering, normalization, and label transfer
