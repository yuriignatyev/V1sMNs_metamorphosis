# convert_cellbender_h5_to_seurat.R

# Increase object size limit for large matrices
options(future.globals.maxSize = 4096 * 1024^2)

# libraries
library(rhdf5)
library(Matrix)
library(Seurat)

# ----- input -----
# replace this path with your own CellBender output file (.h5)
h5file <- "path/to/your/output_filtered.h5"

# ----- function to convert h5 file from cellbender onto seurat data object -----
convert_cellbender_to_seurat <- function(h5file) {
                
     
  barcodes <- h5read(h5file, "/matrix/barcodes")
  data     <- h5read(h5file, "/matrix/data")
  indices  <- h5read(h5file, "/matrix/indices") + 1 # fixes 0-based indexing
  indptr   <- h5read(h5file, "/matrix/indptr")
  shape    <- h5read(h5file, "/matrix/shape")
  genes    <- h5read(h5file, "/matrix/features/name")

  #ensuring unique names
  genes <- make.unique(genes)
  barcodes <- make.unique(barcodes)

  #building sparse matrix
  mat <- sparseMatrix(i = indices, p = indptr, x = data, dims = c(shape[1], shape[2]))  
  rownames(mat) <- genes
  colnames(mat) <- barcodes
                           
  seu <- CreateSeuratObject(mat)
  return(seu)
}            

# ----- creating seurat object -----
seu <- convert_cellbender_to_seurat(h5file)


# saveRDS(seu, file = "your_seurat_object.rds")
