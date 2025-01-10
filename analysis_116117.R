library(Seurat)
library(ggplot2)

#Adult 116/117 (192(old ident))
# Load data
TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
Adult <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/Adult.rds")
# Set default assay and identities
DefaultAssay(Adult) <- "RNA"
Idents(Adult) <- Adult$MultiomeNN
# Check the original distribution of clusters
table(Adult$MultiomeNN)
# Subset cells labeled as 192 within clusters 116 and 117
sub_116 <- subset(Adult, idents = 116)
Idents(sub_116) <- sub_116$FinalIdents
cells_116 <- WhichCells(sub_116, idents = 192)
sub_117 <- subset(Adult, idents = 117)
Idents(sub_117) <- sub_117$FinalIdents
cells_117 <- WhichCells(sub_117, idents = 192)
# Merge these subsets into a new cluster called "0"
cells_to_cluster_0 <- c(cells_116, cells_117)
# Ensure MultiomeNN is a factor and update
Adult$MultiomeNN <- as.character(Adult$MultiomeNN) # Convert to character to allow new cluster addition
Adult$MultiomeNN[cells_to_cluster_0] <- "0"
Adult$MultiomeNN <- factor(Adult$MultiomeNN, levels = c("0", as.character(1:259)))
# Check the new distribution of clusters
table(Adult$MultiomeNN)


#P50 116/117 cells (0 old ident)
P50 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P50.rds")
DefaultAssay(P50) <- "RNA"
table(P50@active.ident)
Idents(P50) <- P50$MultiomeNN
# Display the new distribution of clusters
table(P50$MultiomeNN)
# Subset cells labeled as 192 within clusters 116 and 117
sub_116_P50 <- subset(P50, idents = 116)
Idents(sub_116_P50) <- sub_116_P50$FinalIdents.T45separate
cells_116_P50 <- WhichCells(sub_116_P50, idents = 0)
sub_117_P50 <- subset(P50, idents = 117)
Idents(sub_117_P50) <- sub_117_P50$FinalIdents.T45separate
cells_117_P50 <- WhichCells(sub_117_P50, idents = 0)
# Merge these subsets into a new cluster called "0"
cells_to_cluster_0_P50 <- c(cells_116_P50, cells_117_P50)
# Ensure MultiomeNN is a factor and update
P50$MultiomeNN <- as.character(P50$MultiomeNN) # Convert to character to allow new cluster addition
P50$MultiomeNN[cells_to_cluster_0_P50] <- "0"
P50$MultiomeNN <- factor(P50$MultiomeNN, levels = c("0", as.character(1:259)))
# Display the new distribution of clusters
table(P50$MultiomeNN)


#P30 116/117 ( 0 old identities)

P30 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P30.rds")
# Set default assay and identities
DefaultAssay(P30) <- "RNA"
table(P30@active.ident)  # Display the current cluster distribution
Idents(P30) <- P30$MultiomeNN
# Display the current distribution of clusters
table(P30$MultiomeNN)
# Subset cells labeled as 0 within clusters 116 and 117
sub_116_P30 <- subset(P30, idents = 116)
Idents(sub_116_P30) <- sub_116_P30$FinalIdents
cells_116_P30 <- WhichCells(sub_116_P30, idents = 0)
sub_117_P30 <- subset(P30, idents = 117)
Idents(sub_117_P30) <- sub_117_P30$FinalIdents
cells_117_P30 <- WhichCells(sub_117_P30, idents = 0)
# Merge these subsets into a new cluster called "0"
cells_to_cluster_0_P30 <- c(cells_116_P30, cells_117_P30)
# Ensure MultiomeNN is a factor and update
P30$MultiomeNN <- as.character(P30$MultiomeNN)  # Convert to character to allow new cluster addition
P30$MultiomeNN[cells_to_cluster_0_P30] <- "0"
# Reorder the 'MultiomeNN' factor to ensure it's in the correct order
P30$MultiomeNN <- factor(P30$MultiomeNN, levels = c("0", as.character(1:259)))
# Display the new distribution of clusters after reordering
table(P30$MultiomeNN)