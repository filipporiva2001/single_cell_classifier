library(Seurat)
library(ggplot2)

TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
#P50
P50 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P50.rds")
DefaultAssay(P50) <- "RNA"
table(P50@active.ident)
Idents(P50) <- P50$MultiomeNN
#Cluster 209 and 210 (different due to technical issues)
sub <- subset(P50, idents = 22)
Idents(sub) <- sub$FinalIdents.T45separate
table(Idents(sub))
markers.c22 <- FindMarkers(sub, ident.1 = 209,ident.2 = 210, assay = "RNA")
markers.c22 <- markers.c22[markers.c22$p_val_adj<0.05,]
markers.c22.tf <- markers.c22[rownames(markers.c22)%in%TFs,]

sub <- subset(P50, subset = FinalIdents.T45separate == 209)
table(Idents(sub))
markers.209 <- FindMarkers(sub, ident.1 = 22,ident.2 = 2, assay = "RNA")
markers.209 <- markers.209[markers.209$p_val_adj<0.05,]
markers.209.tf <- markers.209[rownames(markers.209)%in%TFs,]

#Cluster 85 check 200 and 201
sub <- subset(P50, idents = 85)
Idents(sub) <- sub$FinalIdents.T45separate
table(Idents(sub))
markers.c85 <- FindMarkers(sub, ident.1 = 200,ident.2 = 201, assay = "RNA")
markers.c85 <- markers.c85[markers.c85$p_val_adj<0.05,]
markers.c85.tf <- markers.c85[rownames(markers.c85)%in%TFs,]


sub <- subset(P50, subset = FinalIdents.T45separate == 200)
table(Idents(sub))
markers.200 <- FindMarkers(sub, ident.1 = 84,ident.2 = 85, assay = "RNA")
markers.200 <- markers.200[markers.200$p_val_adj<0.05,]
markers.200.tf <- markers.200[rownames(markers.200)%in%TFs,]

sub <- subset(P50, subset = FinalIdents.T45separate == 201)
table(Idents(sub))
markers.201 <- FindMarkers(sub, ident.1 = 84,ident.2 = 85, assay = "RNA")
markers.201 <- markers.201[markers.201$p_val_adj<0.05,]
markers.201.tf <- markers.201[rownames(markers.201)%in%TFs,]

#Cluster 237 check 232

sub <- subset(P50, subset = FinalIdents.T45separate == 232)
table(Idents(sub))
x <- WhichCells(sub, idents = 237, invert = TRUE) #collect all cell that are not 237 to x parameter 
sub <- SetIdent(sub, cells = x, value = 0) #assign to x parameter the value 0
table(Idents(sub))
markers.232 <- FindMarkers(sub, ident.1 = 237,ident.2 = 0, assay = "RNA")
markers.232 <- markers.232[markers.232$p_val_adj<0.05,]
markers.232.tf <- markers.232[rownames(markers.232)%in%TFs,]

#Adult 
library(Seurat)
library(ggplot2)

TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
Adult <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/Adult.rds")
DefaultAssay(Adult) <- "RNA"
table(Adult@active.ident)
Idents(Adult) <- Adult$MultiomeNN

#61 should be classified as 136 not 61, I compared 61 with 60 in 19 they appeared to be different
sub <- subset(Adult, subset = FinalIdents == 19)
table(Idents(sub))
markers.19 <- FindMarkers(sub, ident.1 = 61, ident.2 = 60, assay = "RNA")
markers.19 <- markers.19[markers.19$p_val_adj < 0.05,]
markers.19.tf <- markers.19[rownames(markers.19) %in% TFs,]

#85 Compare 192 with rest of the cells (85 should be classified as 201 not 192)
# Subset Adult object based on FinalIdents
sub <- subset(Adult, idents = 85)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 192, invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 0)
table(Idents(sub))
markers.c85 <- FindMarkers(sub, ident.1 = 192, ident.2 = 0, assay = "RNA")
markers.c85 <- markers.c85[markers.c85$p_val_adj < 0.05,]
markers.c85.tf <- markers.c85[rownames(markers.c85) %in% TFs,]

table(Adult@active.ident)
Idents(Adult) <- Adult$MultiomeNN
markers.c201 <- FindMarkers(Adult, ident.1 = 84, ident.2 = 85, assay = "RNA")
markers.c201 <- markers.c201[markers.c201$p_val_adj < 0.05,]
markers.c201.tf <- markers.c201[rownames(markers.c201) %in% TFs,]
View(markers.c201)
View(markers.c201.tf)
Idents(P50) <- P50$MultiomeNN
markers.P50201 <- FindMarkers(P50, ident.1 = 84, ident.2 = 85, assay = "RNA")
markers.P50201 <- markers.P50201[markers.P50201$p_val_adj < 0.05,]
markers.P50201.tf <- markers.P50201[rownames(markers.P50201) %in% TFs,]
View(markers.P50201)
View(markers.P50201.tf)


#96 Compare 28 with 85 cells 
sub <- subset(Adult, idents = 96)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c96 <- FindMarkers(sub, ident.1 = 28, ident.2 = 85, assay = "RNA")
markers.c96 <- markers.c96[markers.c96$p_val_adj < 0.05,]
markers.c96.tf <- markers.c96[rownames(markers.c96) %in% TFs,]

#116 compare 192 with other cells 
sub <- subset(Adult, idents = 116)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c116 <- FindMarkers(sub, ident.1 = 205, ident.2 = 186, assay = "RNA")
markers.c116 <- markers.c116[markers.c116$p_val_adj < 0.05,]
markers.c116.tf <- markers.c116[rownames(markers.c116) %in% TFs,]

sub <- subset(Adult, idents = 116)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 192, invert = TRUE) #collect all cells that are not 192 to x parameter 
sub <- SetIdent(sub, cells = x, value = 0) #assign to x parameter the value 0
table(Idents(sub))
markers.116 <- FindMarkers(sub, ident.1 = 192, ident.2 = 0, assay = "RNA")
markers.116 <- markers.116[markers.116$p_val_adj < 0.05, ]
markers.116.tf <- markers.116[rownames(markers.116) %in% TFs, ]


#117 compare 192 with all other cells 
sub <- subset(Adult, idents = 117)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c117 <- FindMarkers(sub, ident.1 = 192, ident.2 = 185, assay = "RNA")
markers.c117 <- markers.c117[markers.c117$p_val_adj < 0.05,]
markers.c117.tf <- markers.c117[rownames(markers.c117) %in% TFs,]



sub <- subset(Adult, subset = FinalIdents == 192)
table(Idents(sub))
markers.192 <- FindMarkers(sub, ident.1 = 116, ident.2 = 117, assay = "RNA")
markers.192 <- markers.192[markers.192$p_val_adj < 0.05,]
markers.192.tf <- markers.192[rownames(markers.192) %in% TFs,]


#182 Compare 91 to all other cells 
sub <- subset(Adult, idents = 182)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 91, invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 0)
table(Idents(sub))
markers.c182 <- FindMarkers(sub, ident.1 = 91, ident.2 = 0, assay = "RNA")
markers.c182 <- markers.c182[markers.c182$p_val_adj < 0.05,]
markers.c182.tf <- markers.c182[rownames(markers.c182) %in% TFs,]

#202  Compare 231 with 232 
sub <- subset(Adult, idents = 202)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c202 <- FindMarkers(sub, ident.1 = 231, ident.2 = 232, assay = "RNA")
markers.c202 <- markers.c202[markers.c202$p_val_adj < 0.05,]
markers.c202.tf <- markers.c202[rownames(markers.c202) %in% TFs,]

#237 compare 231 and 232 at 202 

sub <- subset(Adult, idents = 202)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c202 <- FindMarkers(sub, ident.1 = 231, ident.2 = 232, assay = "RNA")
markers.c202 <- markers.c202[markers.c202$p_val_adj < 0.05,]
markers.c202.tf <- markers.c202[rownames(markers.c202) %in% TFs,]

sub <- subset(P50, idents = 202)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c202 <- FindMarkers(P50, ident.1 = 231, ident.2 = 232, assay = "RNA")
markers.c202 <- markers.c202[markers.c202$p_val_adj < 0.05,]
markers.c202.tf <- markers.c202[rownames(markers.c202) %in% TFs,]
#comparison 231/232 at Adult stage 
TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
Adult <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/Adult.rds")
DefaultAssay(Adult) <- "RNA"
table(Adult@active.ident)
Idents(Adult) <- Adult$MultiomeNN
markers.c202 <- FindMarkers(Adult, ident.1 = 231, ident.2 = 232, assay = "RNA")
markers.c202 <- markers.c202[markers.c202$p_val_adj < 0.05,]
markers.c202.tf <- markers.c202[rownames(markers.c202) %in% TFs,]
View(markers.c202)
View(markers.c202.tf)

#Comparison 231/232 at P50 stage 
P50 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P50.rds")
DefaultAssay(P50) <- "RNA"
table(P50@active.ident)
Idents(P50) <- P50$MultiomeNN
markers.P50202 <- FindMarkers(P50, ident.1 = 231, ident.2 = 232, assay = "RNA")
markers.P50202 <- markers.P50202[markers.P50202$p_val_adj < 0.05,]
markers.P50202.tf <- markers.P50202[rownames(markers.P50202) %in% TFs,]
View(markers.P50202)
View(markers.P50202.tf)

#Comparison 202 and 237 at Adult and At P50 
# Comparison 237/202 at Adult stage
TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
Adult <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/Adult.rds")
DefaultAssay(Adult) <- "RNA"
table(Adult@active.ident)
Idents(Adult) <- Adult$MultiomeNN
markers.c237 <- FindMarkers(Adult, ident.1 = 237, ident.2 = 202, assay = "RNA")
markers.c237 <- markers.c237[markers.c237$p_val_adj < 0.05,]
markers.c237.tf <- markers.c237[rownames(markers.c237) %in% TFs,]
View(markers.c237)
View(markers.c237.tf)

# Comparison 237/202 at P50 stage
P50 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P50.rds")
DefaultAssay(P50) <- "RNA"
table(P50@active.ident)
Idents(P50) <- P50$MultiomeNN
markers.P50237 <- FindMarkers(P50, ident.1 = 237, ident.2 = 202, assay = "RNA")
markers.P50237 <- markers.P50237[markers.P50237$p_val_adj < 0.05,]
markers.P50237.tf <- markers.P50237[rownames(markers.P50237) %in% TFs,]
View(markers.P50237)
View(markers.P50237.tf)


#240 compare 100 and 102
sub <- subset(Adult, idents = 240)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c240 <- FindMarkers(sub, ident.1 = 100, ident.2 = 102, assay = "RNA")
markers.c240 <- markers.c240[markers.c240$p_val_adj < 0.05,]
markers.c240.tf <- markers.c240[rownames(markers.c240) %in% TFs,]

sub <- subset(Adult, subset = FinalIdents == 102)
table(Idents(sub))
markers.102 <- FindMarkers(sub, ident.1 = 240, ident.2 = 239, assay = "RNA")
markers.102 <- markers.102[markers.102$p_val_adj < 0.05,]
markers.102.tf <- markers.102[rownames(markers.102) %in% TFs,]

sub <- subset(Adult, subset = FinalIdents == 100)
table(Idents(sub))
markers.100 <- FindMarkers(sub, ident.1 = 240, ident.2 = 239, assay = "RNA")
markers.100 <- markers.100[markers.100$p_val_adj < 0.05,]
markers.100.tf <- markers.100[rownames(markers.100) %in% TFs,]

#Cluster 240 
#FIrst compare o100 with o102 at P50 and adult
#Second compare n240 with n141 at P50 and Adult 
markers.c240 <- FindMarkers(Adult, ident.1 = 100, ident.2 = 102, assay = "RNA")
markers.c240 <- markers.c240[markers.c240$p_val_adj < 0.05,]
markers.c240.tf <- markers.c240[rownames(markers.c240) %in% TFs,]
View(markers.c240)
View(markers.c240.tf)
markers.P50240 <- FindMarkers(P50, ident.1 = 100, ident.2 = 102, assay = "RNA")
markers.P50240 <- markers.P50240[markers.P50240$p_val_adj < 0.05,]
markers.P50240.tf <- markers.P50240[rownames(markers.P50240) %in% TFs,]
View(markers.P50240)
View(markers.P50240.tf)

markers.c102 <- FindMarkers(Adult, ident.1 = 100, ident.2 = 102, assay = "RNA")
markers.c102 <- markers.c102[markers.c102$p_val_adj < 0.05,]
markers.c102.tf <- markers.c102[rownames(markers.c102) %in% TFs,]
View(markers.c102)
View(markers.c102.tf)
Idents(P50) <- P50$MultiomeNN
markers.P50102 <- FindMarkers(P50, ident.1 = 100, ident.2 = 102, assay = "RNA")
markers.P50102 <- markers.P50102[markers.P50102$p_val_adj < 0.05,]
markers.P50102.tf <- markers.P50102[rownames(markers.P50102) %in% TFs,]
View(markers.P50102)
View(markers.P50102.tf)


#P30
library(Seurat)
library(ggplot2)

TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
P30 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P30.rds")
DefaultAssay(P30) <- "RNA"
table(P30@active.ident)
Idents(P30) <- P30$MultiomeNN


# 2 compare 0 and 192 
sub <- subset(P30, idents = 2)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c2 <- FindMarkers(sub, ident.1 = 0, ident.2 = 192, assay = "RNA")
markers.c2 <- markers.c2[markers.c2$p_val_adj < 0.05,]
markers.c2.tf <- markers.c2[rownames(markers.c2) %in% TFs,]
#compare cells 192 and 0 to all otehyr cell in the cluster
sub <- subset(P30, idents = 2)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = c(192, 0), invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 1)  # Assign identity 1 to all cells in x
sub <- SetIdent(sub, cells = WhichCells(sub, idents = 192), value = 0)  # Set identity 0 for cells with idents = 192
table(Idents(sub))
markers.2 <- FindMarkers(sub, ident.1 = 0, ident.2 = 1, assay = "RNA")
markers.2 <- markers.2[markers.2$p_val_adj < 0.05,]
markers.2.tf <- markers.2[rownames(markers.2) %in% TFs,]


#5 Compare 219 with 0  in 5 and then compare 219 in 5 and 6 
sub <- subset(P30, idents = 5)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c5 <- FindMarkers(sub, ident.1 = 0, ident.2 = 219, assay = "RNA")
markers.c5 <- markers.c5[markers.c5$p_val_adj < 0.05,]
markers.c5.tf <- markers.c5[rownames(markers.c5) %in% TFs,]

sub <- subset(P30, subset = FinalIdents == 219)
table(Idents(sub))
markers.219 <- FindMarkers(sub, ident.1 = 5, ident.2 = 6, assay = "RNA")
markers.219 <- markers.219[markers.219$p_val_adj < 0.05,]
markers.219.tf <- markers.219[rownames(markers.219) %in% TFs,]


#46 compare cell 0 in 46 with cells 58 in 66
#Comparison n46 with n66 containing 
# Load transcription factors and P30 dataset
TFs <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/FlyBaseTFsConverted.rds")
P30 <- readRDS("/n/sci/SCI-004375-NYUDATA/Filippo/P30.rds")
DefaultAssay(P30) <- "RNA"
# Set identifiers and perform marker comparison for P30 stage
table(P30@active.ident)
Idents(P30) <- P30$MultiomeNN
markers.c58 <- FindMarkers(P30, ident.1 = 46, ident.2 = 66, assay = "RNA")
markers.c58 <- markers.c58[markers.c58$p_val_adj < 0.05,]
markers.c58.tf <- markers.c58[rownames(markers.c58) %in% TFs,]
View(markers.c58)
View(markers.c58.tf)

markers.P5058 <- FindMarkers(P50, ident.1 = 46, ident.2 = 66, assay = "RNA")
markers.P5058 <- markers.P5058[markers.P5058$p_val_adj < 0.05,]
markers.P5058.tf <- markers.P5058[rownames(markers.P5058) %in% TFs,]
View(markers.P5058)
View(markers.P5058.tf)


#96 Compare 28 to all other cells 
sub <- subset(P30, idents = 96)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 28, invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 0)
table(Idents(sub))
markers.c96 <- FindMarkers(sub, ident.1 = 28, ident.2 = 0, assay = "RNA")
markers.c96 <- markers.c96[markers.c96$p_val_adj < 0.05,]
markers.c96.tf <- markers.c96[rownames(markers.c96) %in% TFs,]

#116 Compare 0 with other cells 
sub <- subset(P30, idents = 116)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 0, invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 1)
table(Idents(sub))
markers.c116 <- FindMarkers(sub, ident.1 = 0, ident.2 = 1, assay = "RNA")
markers.c116 <- markers.c116[markers.c116$p_val_adj < 0.05,]
markers.c116.tf <- markers.c116[rownames(markers.c116) %in% TFs,]

sub <- subset(P30, idents = 116)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c116 <- FindMarkers(sub, ident.1 = 0, ident.2 = 187, assay = "RNA")
markers.c116 <- markers.c116[markers.c116$p_val_adj < 0.05,]
markers.c116.tf <- markers.c116[rownames(markers.c116) %in% TFs,]

#117 Compare 0 with other cells 

sub <- subset(P30, idents = 117)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
x <- WhichCells(sub, idents = 0, invert = TRUE)
sub <- SetIdent(sub, cells = x, value = 1)
table(Idents(sub))
markers.c117 <- FindMarkers(sub, ident.1 = 0, ident.2 = 1, assay = "RNA")
markers.c117 <- markers.c117[markers.c117$p_val_adj < 0.05,]
markers.c117.tf <- markers.c117[rownames(markers.c117) %in% TFs,]

sub <- subset(P30, idents = 117)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c117 <- FindMarkers(sub, ident.1 = 0, ident.2 = 196, assay = "RNA")
markers.c117 <- markers.c117[markers.c117$p_val_adj < 0.05,]
markers.c117.tf <- markers.c117[rownames(markers.c117) %in% TFs,]

#Compare 117 and 116 in 0 
sub <- subset(P30, subset = FinalIdents == 0)
table(Idents(sub))
markers.0 <- FindMarkers(sub, ident.1 = 116, ident.2 = 117, assay = "RNA")
markers.0 <- markers.0[markers.0$p_val_adj < 0.05,]
markers.0.tf <- markers.0[rownames(markers.0) %in% TFs,]

#125 compare 0 and 83 cells 
sub <- subset(P30, idents = 125)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c125 <- FindMarkers(sub, ident.1 = 0, ident.2 = 83, assay = "RNA")
markers.c125 <- markers.c125[markers.c125$p_val_adj < 0.05,]
markers.c125.tf <- markers.c125[rownames(markers.c125) %in% TFs,]

#Comparte n125 and n90
Idents(P30) <- P30$MultiomeNN
markers.c83 <- FindMarkers(P30, ident.1 = 125, ident.2 = 90, assay = "RNA")
markers.c83 <- markers.c83[markers.c83$p_val_adj < 0.05,]
markers.c83.tf <- markers.c83[rownames(markers.c83) %in% TFs,]
View(markers.c83)
View(markers.c83.tf)

Idents(P50) <- P50$MultiomeNN
markers.P5083 <- FindMarkers(P50, ident.1 = 125, ident.2 = 90, assay = "RNA")
markers.P5083 <- markers.P5083[markers.P5083$p_val_adj < 0.05,]
markers.P5083.tf <- markers.P5083[rownames(markers.P5083) %in% TFs,]
View(markers.P5083)
View(markers.P5083.tf)

#240/241 compare if 0 cells are different between the two clusters 
sub <- subset(P30, subset = FinalIdents == 0)
table(Idents(sub))
markers.0 <- FindMarkers(P30, ident.1 = 240, ident.2 = 241, assay = "RNA")
markers.0 <- markers.0[markers.0$p_val_adj < 0.05,]
markers.0.tf <- markers.0[rownames(markers.0) %in% TFs,]
#Cluster 240 and 241 
# Set identifiers and perform marker comparison for P30 stage
Idents(P30) <- P30$MultiomeNN
markers.c0 <- FindMarkers(P30, ident.1 = 240, ident.2 = 241, assay = "RNA")
markers.c0 <- markers.c0[markers.c0$p_val_adj < 0.05,]
markers.c0.tf <- markers.c0[rownames(markers.c0) %in% TFs,]
View(markers.c0)
View(markers.c0.tf)

#258 compare 67 with 0 
sub <- subset(P30, idents = 258)
Idents(sub) <- sub$FinalIdents
table(Idents(sub))
markers.c258 <- FindMarkers(sub, ident.1 = 0, ident.2 = 67, assay = "RNA")
markers.c258 <- markers.c258[markers.c258$p_val_adj < 0.05,]
markers.c258.tf <- markers.c258[rownames(markers.c258) %in% TFs,]


