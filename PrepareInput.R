library(Seurat)

# Build a Seurat object with your dataset.

# Log-normalize to 10,000
Adult.male <- NormalizeData(Adult.male, normalization.method = "LogNormalize", scale.factor = 10000)

# Load the set of markers for the stage closest to the target dataset.
markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNAdult.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNP70.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNP50.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNP40.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNP30.rds")
# markers <- readRDS("C:/YourDirectory/AppendixNN/MarkersNNP15.rds")

# Convert the marker names to match your dataset first!

# Subset for the marker expression and save as a full matrix.
Adult.male.markers <- as.matrix(Adult.male@assays$RNA@data[markers,])
saveRDS(Adult.male.markers, file = "Example.rds")
