setwd("/n/sci/SCI-004375-NYUDATA/Filippo")

library(Seurat)
library(dplyr)


# Load annotation data
annotation <- read.csv(file = "AnnotationJul24.csv")
annotated <- annotation[which(nchar(annotation$annotation)!=0),]
celltypes <- annotated$annotation
names(celltypes) <- annotated$cluster

#EarlyP24
Zipursky.Early <- readRDS("Zipursky_Early.rds")
Preds.Zipursky.Early <- read.table(file = "EarlyZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for Early
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.Early) <- colnames(Zipursky.Early@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.Early) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.Early <- AddMetaData(Zipursky.Early, metadata = Preds.Zipursky.Early)
# Run the summary() command right after adding the metadata
summary(Zipursky.Early@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.Early <- SetIdent(Zipursky.Early, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.Early <- RenameIdents(Zipursky.Early, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.Early$MultiomeAnnotated <- Idents(Zipursky.Early)
# Save the modified Seurat object
saveRDS(Zipursky.Early, file = "Zipursky_Earlyupdate.rds")

Zipursky.Early <- FindVariableFeatures(Zipursky.Early)
Zipursky.Early <- ScaleData(Zipursky.Early)
Zipursky.Early <- RunPCA(Zipursky.Early, npcs = 100)
Zipursky.Early <- RunTSNE(Zipursky.Early, dims = 1:100)


DimPlot(Zipursky.Early, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.Early, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.Early, group.by = "set")


#P24P24
Zipursky.P24 <- readRDS("Zipursky_P24.rds")
Preds.Zipursky.P24 <- read.table(file = "P24ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P24
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P24) <- colnames(Zipursky.P24@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P24) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P24 <- AddMetaData(Zipursky.P24, metadata = Preds.Zipursky.P24)
# Run the summary() command right after adding the metadata
summary(Zipursky.P24@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P24 <- SetIdent(Zipursky.P24, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P24 <- RenameIdents(Zipursky.P24, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P24$MultiomeAnnotated <- Idents(Zipursky.P24)
# Save the modified Seurat object
saveRDS(Zipursky.P24, file = "Zipursky_P24update.rds")

Zipursky.P24 <- FindVariableFeatures(Zipursky.P24)
Zipursky.P24 <- ScaleData(Zipursky.P24)
Zipursky.P24 <- RunPCA(Zipursky.P24, npcs = 100)
Zipursky.P24 <- RunTSNE(Zipursky.P24, dims = 1:100)

#DimPlot(Zipursky.P24, label=TRUE) +NoLegend()
#FeaturePlot(Zipursky.P24, features = 'MultiomeNN_Conf')
#DimPlot(Zipursky.P24, group.by = "set") 

#P36P24
Zipursky.P36 <- readRDS("Zipursky_P36.rds")
Preds.Zipursky.P36 <- read.table(file = "P36P24ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P36
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P36) <- colnames(Zipursky.P36@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P36) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P36 <- AddMetaData(Zipursky.P36, metadata = Preds.Zipursky.P36)
# Run the summary() command right after adding the metadata
summary(Zipursky.P36@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P36 <- SetIdent(Zipursky.P36, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P36 <- RenameIdents(Zipursky.P36, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P36$MultiomeAnnotated <- Idents(Zipursky.P36)
# Save the modified Seurat object
saveRDS(Zipursky.P36, file = "Zipursky_P36P24update.rds")

Zipursky.P36 <- FindVariableFeatures(Zipursky.P36)
Zipursky.P36 <- ScaleData(Zipursky.P36)
Zipursky.P36 <- RunPCA(Zipursky.P36, npcs = 100)
Zipursky.P36 <- RunTSNE(Zipursky.P36, dims = 1:100)

DimPlot(Zipursky.P36, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P36, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P36, group.by = "set") 



#P36P48 
Zipursky.P36 <- readRDS("Zipursky_P36.rds")
Preds.Zipursky.P36 <- read.table(file = "P36P48ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P36
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P36) <- colnames(Zipursky.P36@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P36) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P36 <- AddMetaData(Zipursky.P36, metadata = Preds.Zipursky.P36)
# Run the summary() command right after adding the metadata
summary(Zipursky.P36@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P36 <- SetIdent(Zipursky.P36, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P36 <- RenameIdents(Zipursky.P36, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P36$MultiomeAnnotated <- Idents(Zipursky.P36)
# Save the modified Seurat object
saveRDS(Zipursky.P36, file = "Zipursky_P36P48update.rds")

Zipursky.P36 <- FindVariableFeatures(Zipursky.P36)
Zipursky.P36 <- ScaleData(Zipursky.P36)
Zipursky.P36 <- RunPCA(Zipursky.P36, npcs = 100)
Zipursky.P36 <- RunTSNE(Zipursky.P36, dims = 1:100)

DimPlot(Zipursky.P36, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P36, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P36, group.by = "set")


#P48P48
#P48
Zipursky.P48 <- readRDS("Zipursky_P48.rds")
Preds.Zipursky.P48 <- read.table(file = "P48ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P48
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P48) <- colnames(Zipursky.P48@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P48) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P48 <- AddMetaData(Zipursky.P48, metadata = Preds.Zipursky.P48)
# Run the summary() command right after adding the metadata
summary(Zipursky.P48@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P48 <- SetIdent(Zipursky.P48, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P48 <- RenameIdents(Zipursky.P48, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P48$MultiomeAnnotated <- Idents(Zipursky.P48)
# Save the modified Seurat object
saveRDS(Zipursky.P48, file = "Zipursky_P48update.rds")

Zipursky.P48 <- FindVariableFeatures(Zipursky.P48)
Zipursky.P48 <- ScaleData(Zipursky.P48)
Zipursky.P48 <- RunPCA(Zipursky.P48, npcs = 100)
Zipursky.P48 <- RunTSNE(Zipursky.P48, dims = 1:100)

#DimPlot(Zipursky.P48, label=TRUE) +NoLegend()
#FeaturePlot(Zipursky.P48, features = 'MultiomeNN_Conf')
#DimPlot(Zipursky.P48, group.by = "set")

#P60P48
Zipursky.P60 <- readRDS("Zipursky_P60.rds")
Preds.Zipursky.P60 <- read.table(file = "P60P48ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P60
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P60) <- colnames(Zipursky.P60@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P60) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P60 <- AddMetaData(Zipursky.P60, metadata = Preds.Zipursky.P60)
# Run the summary() command right after adding the metadata
summary(Zipursky.P60@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P60 <- SetIdent(Zipursky.P60, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P60 <- RenameIdents(Zipursky.P60, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P60$MultiomeAnnotated <- Idents(Zipursky.P60)
# Save the modified Seurat object
saveRDS(Zipursky.P60, file = "Zipursky_P60P48update.rds")

Zipursky.P60 <- FindVariableFeatures(Zipursky.P60)
Zipursky.P60 <- ScaleData(Zipursky.P60)
Zipursky.P60 <- RunPCA(Zipursky.P60, npcs = 100)
Zipursky.P60 <- RunTSNE(Zipursky.P60, dims = 1:100)

DimPlot(Zipursky.P60, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P60, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P60, group.by = "set")

#P60Adult
Zipursky.P60 <- readRDS("Zipursky_P60.rds")
Preds.Zipursky.P60 <- read.table(file = "P60AdultZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P60
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P60) <- colnames(Zipursky.P60@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P60) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P60 <- AddMetaData(Zipursky.P60, metadata = Preds.Zipursky.P60)
# Run the summary() command right after adding the metadata
summary(Zipursky.P60@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P60 <- SetIdent(Zipursky.P60, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P60 <- RenameIdents(Zipursky.P60, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P60$MultiomeAnnotated <- Idents(Zipursky.P60)
# Save the modified Seurat object
saveRDS(Zipursky.P60, file = "Zipursky_P60Adultupdate.rds")

Zipursky.P60 <- FindVariableFeatures(Zipursky.P60)
Zipursky.P60 <- ScaleData(Zipursky.P60)
Zipursky.P60 <- RunPCA(Zipursky.P60, npcs = 100)
Zipursky.P60 <- RunTSNE(Zipursky.P60, dims = 1:100)

DimPlot(Zipursky.P60, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P60, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P60, group.by = "set")


#P72P48
Zipursky.P72 <- readRDS("Zipursky_P72.rds")
Preds.Zipursky.P72 <- read.table(file = "P72P48ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P72
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P72) <- colnames(Zipursky.P72@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P72) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P72 <- AddMetaData(Zipursky.P72, metadata = Preds.Zipursky.P72)
# Run the summary() command right after adding the metadata
summary(Zipursky.P72@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P72 <- SetIdent(Zipursky.P72, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P72 <- RenameIdents(Zipursky.P72, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P72$MultiomeAnnotated <- Idents(Zipursky.P72)
# Save the modified Seurat object
saveRDS(Zipursky.P72, file = "Zipursky_P72P48update.rds")

Zipursky.P72 <- FindVariableFeatures(Zipursky.P72)
Zipursky.P72 <- ScaleData(Zipursky.P72)
Zipursky.P72 <- RunPCA(Zipursky.P72, npcs = 100)
Zipursky.P72 <- RunTSNE(Zipursky.P72, dims = 1:100)

DimPlot(Zipursky.P72, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P72, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P72, group.by = "set")

#P72Adult
Zipursky.P72 <- readRDS("Zipursky_P72.rds")
Preds.Zipursky.P72 <- read.table(file = "P72AdultZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P72
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P72) <- colnames(Zipursky.P72@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P72) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P72 <- AddMetaData(Zipursky.P72, metadata = Preds.Zipursky.P72)
# Run the summary() command right after adding the metadata
summary(Zipursky.P72@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P72 <- SetIdent(Zipursky.P72, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P72 <- RenameIdents(Zipursky.P72, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P72$MultiomeAnnotated <- Idents(Zipursky.P72)
# Save the modified Seurat object
#saveRDS(Zipursky.P72, file = "Zipursky_P72Adultupdate.rds")

Zipursky.P72 <- FindVariableFeatures(Zipursky.P72)
Zipursky.P72 <- ScaleData(Zipursky.P72)
Zipursky.P72 <- RunPCA(Zipursky.P72, npcs = 100)
Zipursky.P72 <- RunTSNE(Zipursky.P72, dims = 1:100)

#DimPlot(Zipursky.P72, label=TRUE) +NoLegend()
#FeaturePlot(Zipursky.P72, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P72, group.by = "set")

#P84P48
Zipursky.P84 <- readRDS("Zipursky_P84.rds")
Preds.Zipursky.P84 <- read.table(file = "P84P48ZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P84
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P84) <- colnames(Zipursky.P84@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P84) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P84 <- AddMetaData(Zipursky.P84, metadata = Preds.Zipursky.P84)
# Run the summary() command right after adding the metadata
summary(Zipursky.P84@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P84 <- SetIdent(Zipursky.P84, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P84 <- RenameIdents(Zipursky.P84, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P84$MultiomeAnnotated <- Idents(Zipursky.P84)
# Save the modified Seurat object
saveRDS(Zipursky.P84, file = "Zipursky_P84P48update.rds")

Zipursky.P84 <- FindVariableFeatures(Zipursky.P84)
Zipursky.P84 <- ScaleData(Zipursky.P84)
Zipursky.P84 <- RunPCA(Zipursky.P84, npcs = 100)
Zipursky.P84 <- RunTSNE(Zipursky.P84, dims = 1:100)

DimPlot(Zipursky.P84, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P84, features = 'MultiomeNN_Conf')
DimPlot(Zipursky.P84, group.by = "set")

#P84Adult
Zipursky.P84 <- readRDS("Zipursky_P84.rds")
Preds.Zipursky.P84 <- read.table(file = "P84AdultZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P84
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P84) <- colnames(Zipursky.P84@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P84) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P84 <- AddMetaData(Zipursky.P84, metadata = Preds.Zipursky.P84)
# Run the summary() command right after adding the metadata
summary(Zipursky.P84@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P84 <- SetIdent(Zipursky.P84, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P84 <- RenameIdents(Zipursky.P84, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P84$MultiomeAnnotated <- Idents(Zipursky.P84)
# Save the modified Seurat object
#saveRDS(Zipursky.P84, file = "Zipursky_P84Adultupdate.rds")

Zipursky.P84 <- FindVariableFeatures(Zipursky.P84)
Zipursky.P84 <- ScaleData(Zipursky.P84)
Zipursky.P84 <- RunPCA(Zipursky.P84, npcs = 100)
Zipursky.P84 <- RunTSNE(Zipursky.P84, dims = 1:100)

DimPlot(Zipursky.P84, label=TRUE) +NoLegend()
FeaturePlot(Zipursky.P84, features = 'MultiomeNN_Conf')
#DimPlot(Zipursky.P84, group.by = "set")


#P96Adult
Zipursky.P96 <- readRDS("Zipursky_P96.rds")
Preds.Zipursky.P96 <- read.table(file = "P96AdultZipurskyPreds_and_Confidence.txt", sep = ",", header = TRUE)
# Load predictions and confidence data for P96
# Ensure the row names in the predictions match the cell names in the Seurat object
row.names(Preds.Zipursky.P96) <- colnames(Zipursky.P96@assays$RNA@data)
# Rename columns of the predictions file
colnames(Preds.Zipursky.P96) <- c("MultiomeNN", "MultiomeNN_Conf")
# Add metadata (predictions) to Seurat object
Zipursky.P96 <- AddMetaData(Zipursky.P96, metadata = Preds.Zipursky.P96)
# Run the summary() command right after adding the metadata
summary(Zipursky.P96@meta.data$MultiomeNN_Conf)
# Set active identities based on MultiomeNN
Zipursky.P96 <- SetIdent(Zipursky.P96, value = "MultiomeNN")
# Rename identities using cell types from annotation
Zipursky.P96 <- RenameIdents(Zipursky.P96, celltypes)
# Set the renamed identities as a new column in the Seurat object
Zipursky.P96$MultiomeAnnotated <- Idents(Zipursky.P96)
# Save the modified Seurat object
saveRDS(Zipursky.P96, file = "Zipursky_P96update.rds")

Zipursky.P96 <- FindVariableFeatures(Zipursky.P96)
Zipursky.P96 <- ScaleData(Zipursky.P96)
Zipursky.P96 <- RunPCA(Zipursky.P96, npcs = 100)
Zipursky.P96 <- RunTSNE(Zipursky.P96, dims = 1:100)

#DimPlot(Zipursky.P96, label=TRUE) +NoLegend()
#FeaturePlot(Zipursky.P96, features = 'MultiomeNN_Conf')
#DimPlot(Zipursky.P96, group.by = "set")
