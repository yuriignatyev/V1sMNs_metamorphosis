# Label transfer pipeline: annotate spinal neurons using coarse and cardinal class neuronal references

library(Seurat)

# --- Coarse label transfer (e.g. neuron, MN, floor plate, roof plate) ---

anchors <- FindTransferAnchors(
  reference = all,  # coarse reference object
  query = seu,
  normalization.method = "LogNormalize",
  dims = 1:50
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(all),
  prediction.assay = TRUE
)

seu[["predictions"]] <- predictions.assay
seu$prediction.id <- GetTransferPredictions(seu, score.filter = 0.5)
Idents(seu) <- seu$prediction.id

saveRDS(seu, "SeuratObjects/AllCells_after_LT_coarse.rds")

# --- Subset neurons (e.g. Interneurons), reprocess, and cluster ---

ins <- subset(seu, idents = "Interneurons")

ins <- NormalizeData(ins)
ins <- FindVariableFeatures(ins, selection.method = "vst", nfeatures = 7000)
ins <- ScaleData(ins)
ins <- RunPCA(ins, npcs = 50)
ins <- RunUMAP(ins, dims = 1:50)

ins1 <- FindNeighbors(ins, dims = 1:50)
ins1 <- FindClusters(ins1, resolution = 2)

# Remove unwanted clusters (e.g. oligodendrocytes / doublets)
ins1 <- subset(ins1, idents = c(0, 31, 32), invert = TRUE)

# Reprocess after cleanup
ins1 <- NormalizeData(ins1)
ins1 <- FindVariableFeatures(ins1, selection.method = "vst", nfeatures = 7000)
ins1 <- ScaleData(ins1)
ins1 <- RunPCA(ins1, npcs = 50)
ins1 <- RunUMAP(ins1, dims = 1:50)

# --- Fine label transfer onto INs (cardinal classes) ---

ins_ref <- readRDS("References/INs_cardinal_reference.rds")

ins_ref <- NormalizeData(ins_ref)
ins_ref <- FindVariableFeatures(ins_ref, selection.method = "vst", nfeatures = 7000)
ins_ref <- ScaleData(ins_ref)
ins_ref <- RunPCA(ins_ref, npcs = 50)

anchors <- FindTransferAnchors(
  reference = ins_ref,
  query = ins1,
  normalization.method = "LogNormalize",
  dims = 1:50
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(ins_ref),
  prediction.assay = TRUE
)

ins1[["predictions"]] <- predictions.assay
ins1$prediction.id <- GetTransferPredictions(ins1, score.filter = 0.5)
Idents(ins1) <- ins1$prediction.id

saveRDS(ins1, "SeuratObjects/INs_after_LT_fine.rds")

# --- Remove unassigned and finalize INs object ---

ins1 <- subset(ins1, idents = "Unassigned", invert = TRUE)

ins1 <- NormalizeData(ins1)
ins1 <- FindVariableFeatures(ins1, selection.method = "vst", nfeatures = 7000)
ins1 <- ScaleData(ins1)
ins1 <- RunPCA(ins1, npcs = 50)
ins1 <- RunUMAP(ins1, dims = 1:50)

saveRDS(ins1, "SeuratObjects/INs_cardinal_cleaned.rds")

# --- Merge INs and MNs for full neuronal object ---

mns <- subset(seu, idents = "Motor Neurons")
neu <- merge(ins1, mns)

neu <- NormalizeData(neu)
neu <- FindVariableFeatures(neu, selection.method = "vst", nfeatures = 7000)
neu <- ScaleData(neu)
neu <- RunPCA(neu, npcs = 50)
neu <- RunUMAP(neu, dims = 1:50)

# Rename MNs label
neu <- RenameIdents(neu, "Motor Neurons" = "MNs")

# Cardinal gene panel (for visualization)
genes <- c(
  "barhl1.L", "lhx2.S", "lhx9.S", "hmx3.L", "foxd3.L", "isl1.L", "tlx3.L",
  "lmx1b.1.S", "sall3.L", "wt1.L", "dmrt3.L", "evx1.L", "evx2.S",
  "en1.L", "foxp2.L", "lhx3.L", "vsx2.S", "shox2.S",
  "sox14.L", "sox21.L", "tal1.S", "gata3.L", "nkx2-2.L", "sim1.S"
)

# Set identity order
levels(neu) <- rev(c("dI1", "dI2", "dI3", "dI5", "dI4", "dI6", "V0", "V1", "V2a", "V2b", "V3", "MNs"))

# Final plots
DimPlot(neu, label = TRUE, label.size = 5, cols = c(
  "dI1" = "orange", "dI2" = "orangered", "dI3" = "tomato1", "dI4" = "goldenrod3",
  "dI5" = "olivedrab2", "dI6" = "darkolivegreen", "V0" = "lightslateblue",
  "V1" = "maroon1", "V2a" = "violet", "V2b" = "violetred", "V3" = "lightpink1", "MNs" = "green3"
)) +
  DotPlot(neu, features = c(genes, "slit2.S", "slc18a3.S", "gad2.L", "slc17a7.L"), 
          cols = c("white", "black"), dot.scale = 8) +
  RotatedAxis()

# Save final object
saveRDS(neu, "SeuratObjects/Neurons_MNs_INs_merged.rds")
