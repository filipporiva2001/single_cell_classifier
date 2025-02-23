library(Seurat)

# Build a Seurat object with your dataset.

# Log-normalize to 10,000
Adult.male <- NormalizeData(Adult.male, normalization.method = "LogNormalize", scale.factor = 10000)

# Load the set of markers for the stage closest to the target dataset.
markers <- readRDS("C:/YourDirectory/AppendixNN/Adult_RNA-NNmarkers.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/P48_RNA-NNmarkers.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/P24_RNA-NNmarkers.rds")

# Convert the marker names to match your dataset first!

# Subset for the marker expression and save as a full matrix.
Adult.male.markers <- as.matrix(Adult.male@assays$RNA@data[markers,])
saveRDS(Adult.male.markers, file = "Example.rds")
