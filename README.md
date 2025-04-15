This repository contains code used for the analysis of the single-cell data from the manuscript  
**"Emergence of limb movement requires multifold increase in spinal inhibitory cell types."**  

Here, **I focus on identifying spinal neural types across _Xenopus_ metamorphosis**.  
The code is designed to accompany the data from the manuscript and enables reproducibility of key downstream processing steps — including ambient RNA removal, quality filtering, label transfer, and figure generation — rather than serving as a full raw-to-clusters preprocessing pipeline.


### Overview

Four NF stages of _X.laevis_ metamorphosis were analyzed: 38, 47, 54, and 57. The pipeline includes:

1. Ambient RNA removal with CellBender
2. **Loading and Converting** CellBender-cleaned data (`*.h5` files) to Seurat objects
3. **Filtering cells** using Gaussian mixture model-based thresholds on `nFeature_RNA` and `nCount_RNA` 
4. **Label transfer**: on a coarse and neuronal level using previously assigned annotations  
5. **Visualization** of UMAPs, dot plots, and neural type proportions 

### Scripts

- `cellbender.sh`: runs CellBender  on `.h5` raw matrix outputs from cellranger count
- `convert_cellbender_h5_to_seurat.R`: loads CellBender `.h5` output and converts it to a Seurat object
- `filter_bimodal_thresholding.R`: applies `mixtools` to determine cutoff thresholds for filtering cells
- `labelling_coarse_and_neuronal.R`: outlines the logic for identifying neurons, annotating neural subtypes, and performing label transfer onto other stages
- `V1MNs_visualizations.R`: generates plots used  in main and supplementary figures






