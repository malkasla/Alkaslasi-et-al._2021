library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

# Import data
cervical_deep.data <- Read10X(data.dir = "cervical")
lumbar_deep.data <- Read10X(data.dir = "lumbar")
thoracic_deep.data <- Read10X(data.dir = "thoracic")

# Preprocessing
colnames(thoracic_deep.data) <- paste0(colnames(thoracic_deep.data),"-thoracic")
thoracics_deep <-CreateSeuratObject(thoracic_deep.data)
thoracics_deep <- PercentageFeatureSet(object = thoracics_deep, pattern = "^mt-", col.name = "percent.mt")
thoracics_deep <- SCTransform(object = thoracics_deep, vars.to.regress = "percent.mt", verbose = FALSE)
thoracics_deep <- RunPCA(object = thoracics_deep, npcs = 30)
ElbowPlot(thoracics_deep, ndims = 30)
thoracics_deep <- RunUMAP(object = thoracics_deep, dims = 1:30, verbose = FALSE)
thoracics_deep <- FindNeighbors(object = thoracics_deep, dims = 1:30, verbose = FALSE)
thoracics_deep <- FindClusters(object = thoracics_deep, verbose = FALSE)
DimPlot(object = thoracics_deep, label = TRUE) + NoLegend()

colnames(lumbar_deep.data) <- paste0(colnames(lumbar_deep.data),"-lumbar")
lumbars_deep <-CreateSeuratObject(lumbar_deep.data)
lumbars_deep <- PercentageFeatureSet(object = lumbars_deep, pattern = "^mt-", col.name = "percent.mt")
lumbars_deep <- SCTransform(object = lumbars_deep, vars.to.regress = "percent.mt", verbose = FALSE)
lumbars_deep <- RunPCA(object = lumbars_deep, npcs = 30)
ElbowPlot(lumbars_deep, ndims = 30)
lumbars_deep <- RunUMAP(object = lumbars_deep, dims = 1:30, verbose = FALSE)
lumbars_deep <- FindNeighbors(object = lumbars_deep, dims = 1:30, verbose = FALSE)
lumbars_deep <- FindClusters(object = lumbars_deep, verbose = FALSE)
DimPlot(object = lumbars_deep, label = TRUE) + NoLegend()

colnames(cervical_deep.data) <- paste0(colnames(cervical_deep.data),"-cervical")
cervicals_deep <- CreateSeuratObject(cervical_deep.data)
cervicals_deep <- PercentageFeatureSet(object = cervicals_deep, pattern = "^mt-", col.name = "percent.mt")
cervicals_deep <- SCTransform(object = cervicals_deep, vars.to.regress = "percent.mt", verbose = FALSE)
cervicals_deep <- RunPCA(object = cervicals_deep, npcs = 30)
ElbowPlot(cervicals_deep, ndims = 30)
cervicals_deep <- RunUMAP(object = cervicals_deep, dims = 1:30, verbose = FALSE)
cervicals_deep <- FindNeighbors(object = cervicals_deep, dims = 1:30, verbose = FALSE)
cervicals_deep <- FindClusters(object = cervicals_deep, verbose = FALSE)
DimPlot(object = cervicals_deep, label = TRUE) + NoLegend()

thoracics_deep$expt <- "Thoracic"
cervicals_deep$expt <- "Cervical"
lumbars_deep$expt <- "Lumbar"

# Import data from M and F runs
cervical_M.data <- Read10X(data.dir = "cervical M")
lumbar_M.data <- Read10X(data.dir = "lumbar M")
thoracic_M.data <- Read10X(data.dir = "thoracic M")
cervical_F.data <- Read10X(data.dir = "cervical F")
lumbar_F.data <- Read10X(data.dir = "lumbar F")
thoracic_F.data <- Read10X(data.dir = "thoracic F")

# Preprocessing
colnames(thoracic_M.data) <- paste0(colnames(thoracic_M.data),"-thoracic")
thoracics_M <-CreateSeuratObject(thoracic_M.data)
thoracics_M <- PercentageFeatureSet(object = thoracics_M, pattern = "^mt-", col.name = "percent.mt")
thoracics_M <- SCTransform(object = thoracics_M, vars.to.regress = "percent.mt", verbose = FALSE)
thoracics_M <- RunPCA(object = thoracics_M, npcs = 30)
ElbowPlot(thoracics_M, ndims = 30)
thoracics_M <- RunUMAP(object = thoracics_M, dims = 1:30, verbose = FALSE)
thoracics_M <- FindNeighbors(object = thoracics_M, dims = 1:30, verbose = FALSE)
thoracics_M <- FindClusters(object = thoracics_M, verbose = FALSE)
DimPlot(object = thoracics_M, label = TRUE) + NoLegend()

colnames(thoracic_F.data) <- paste0(colnames(thoracic_F.data),"-thoracic")
thoracics_F <-CreateSeuratObject(thoracic_F.data)
thoracics_F <- PercentageFeatureSet(object = thoracics_F, pattern = "^mt-", col.name = "percent.mt")
thoracics_F <- SCTransform(object = thoracics_F, vars.to.regress = "percent.mt", verbose = FALSE)
thoracics_F <- RunPCA(object = thoracics_F, npcs = 30)
ElbowPlot(thoracics_F, ndims = 30)
thoracics_F <- RunUMAP(object = thoracics_F, dims = 1:30, verbose = FALSE)
thoracics_F <- FindNeighbors(object = thoracics_F, dims = 1:30, verbose = FALSE)
thoracics_F <- FindClusters(object = thoracics_F, verbose = FALSE)
DimPlot(object = thoracics_F, label = TRUE) + NoLegend()

colnames(cervical_M.data) <- paste0(colnames(cervical_M.data),"-cervical")
cervicals_M <-CreateSeuratObject(cervical_M.data)
cervicals_M <- PercentageFeatureSet(object = cervicals_M, pattern = "^mt-", col.name = "percent.mt")
cervicals_M <- SCTransform(object = cervicals_M, vars.to.regress = "percent.mt", verbose = FALSE)
cervicals_M <- RunPCA(object = cervicals_M, npcs = 30)
ElbowPlot(cervicals_M, ndims = 30)
cervicals_M <- RunUMAP(object = cervicals_M, dims = 1:30, verbose = FALSE)
cervicals_M <- FindNeighbors(object = cervicals_M, dims = 1:30, verbose = FALSE)
cervicals_M <- FindClusters(object = cervicals_M, verbose = FALSE)
DimPlot(object = cervicals_M, label = TRUE) + NoLegend()

colnames(cervical_F.data) <- paste0(colnames(cervical_F.data),"-cervical")
cervicals_F <-CreateSeuratObject(cervical_F.data)
cervicals_F <- PercentageFeatureSet(object = cervicals_F, pattern = "^mt-", col.name = "percent.mt")
cervicals_F <- SCTransform(object = cervicals_F, vars.to.regress = "percent.mt", verbose = FALSE)
cervicals_F <- RunPCA(object = cervicals_F, npcs = 30)
ElbowPlot(cervicals_F, ndims = 30)
cervicals_F <- RunUMAP(object = cervicals_F, dims = 1:30, verbose = FALSE)
cervicals_F <- FindNeighbors(object = cervicals_F, dims = 1:30, verbose = FALSE)
cervicals_F <- FindClusters(object = cervicals_F, verbose = FALSE)
DimPlot(object = cervicals_F, label = TRUE) + NoLegend()

colnames(lumbar_M.data) <- paste0(colnames(lumbar_M.data),"-lumbar")
lumbars_M <-CreateSeuratObject(lumbar_M.data)
lumbars_M <- PercentageFeatureSet(object = lumbars_M, pattern = "^mt-", col.name = "percent.mt")
lumbars_M <- SCTransform(object = lumbars_M, vars.to.regress = "percent.mt", verbose = FALSE)
lumbars_M <- RunPCA(object = lumbars_M, npcs = 30)
ElbowPlot(lumbars_M, ndims = 30)
lumbars_M <- RunUMAP(object = lumbars_M, dims = 1:30, verbose = FALSE)
lumbars_M <- FindNeighbors(object = lumbars_M, dims = 1:30, verbose = FALSE)
lumbars_M <- FindClusters(object = lumbars_M, verbose = FALSE)
DimPlot(object = lumbars_M, label = TRUE) + NoLegend()

colnames(lumbar_F.data) <- paste0(colnames(lumbar_F.data),"-lumbar")
lumbars_F <-CreateSeuratObject(lumbar_F.data)
lumbars_F <- PercentageFeatureSet(object = lumbars_F, pattern = "^mt-", col.name = "percent.mt")
lumbars_F <- SCTransform(object = lumbars_F, vars.to.regress = "percent.mt", verbose = FALSE)
lumbars_F <- RunPCA(object = lumbars_F, npcs = 30)
ElbowPlot(lumbars_F, ndims = 30)
lumbars_F <- RunUMAP(object = lumbars_F, dims = 1:30, verbose = FALSE)
lumbars_F <- FindNeighbors(object = lumbars_F, dims = 1:30, verbose = FALSE)
lumbars_F <- FindClusters(object = lumbars_F, verbose = FALSE)
DimPlot(object = lumbars_F, label = TRUE) + NoLegend()

thoracics_M$sex <- "M"
cervicals_M$sex <- "M"
lumbars_M$sex <- "M"

thoracics_F$sex <- "F"
cervicals_F$sex <- "F"
lumbars_F$sex <- "F"

thoracics_M$expt <- "Thoracic"
cervicals_M$expt <- "Cervical"
lumbars_M$expt <- "Lumbar"

thoracics_F$expt <- "Thoracic"
cervicals_F$expt <- "Cervical"
lumbars_F$expt <- "Lumbar"

# Integrate all 9 datasets
control.anchors_MFdeep <-FindIntegrationAnchors(object.list = list(thoracics_M, thoracics_F, thoracics_deep, cervicals_M, cervicals_F, cervicals_deep, lumbars_M, lumbars_F, lumbars_deep), dims = 1:30)
motor_MFdeep <- IntegrateData(anchorset = control.anchors_MFdeep, dims = 1:30)
DefaultAssay(object = motor_MFdeep) <- "integrated"

# Process integrated data into object that includes all sequenced cells
motor_MFdeep <- ScaleData(object = motor_MFdeep, verbose = FALSE)
motor_MFdeep <- RunPCA(object = motor_MFdeep, npcs = 30, verbose = FALSE)
ElbowPlot(motor_MFdeep, ndims = 40)
motor_MFdeep <- RunUMAP(object = motor_MFdeep, reduction = "pca", dims = 1:30)
motor_MFdeep <- FindNeighbors(object = motor_MFdeep, reduction = "pca", dims = 1:30)
motor_MFdeep <- FindClusters(motor_MFdeep, resolution = 0.8)

FeaturePlot(motor_MFdeep, c("Chat"), label = T)

# Separate cholinergic clusters by Chat expression
mainlymotor_MFdeep <- subset(motor_MFdeep, ident = c(3, 35, 12, 2, 32, 36, 9, 4, 29, 11, 32, 14, 30, 27, 16, 0, 26, 18, 20, 24, 28, 34, 37), invert = T)
DefaultAssay(object = mainlymotor_MFdeep) <- "integrated"
mainlymotor_MFdeep <- RunPCA(object = mainlymotor_MFdeep,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(mainlymotor_MFdeep,ndims=40)
mainlymotor_MFdeep <- RunUMAP(object = mainlymotor_MFdeep, reduction = "pca", dims = c(1:30))
mainlymotor_MFdeep <- FindNeighbors(object = mainlymotor_MFdeep, reduction = "pca", dims = c(1:30))
mainlymotor_MFdeep <- FindClusters(mainlymotor_MFdeep, resolution = 0.6)
mainlymotor_MFdeep <- RunTSNE(object = mainlymotor_MFdeep, reduction = "pca", dims = c(1:30))
DefaultAssay(object = mainlymotor_MFdeep) <- "SCT"

save(mainlymotor_MFdeep, file = "mainlymotor_MFdeep.robj")

# Add cholinergic neuron subtypes to mainly motor metadata
mainlymotor_MFdeep$cholinergictypes <- Idents(mainlymotor_MFdeep)
current_ids <- c(0:22)
mainlymotor_MFdeep@meta.data$cholinergictypes <- plyr::mapvalues(x = mainlymotor_MFdeep@meta.data$cholinergictypes, from = current_ids, to = c("Skeletal Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Skeletal Motor Neurons",
                                                                                                                                               "Skeletal Motor Neurons",
                                                                                                                                               "Skeletal Motor Neurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Visceral Motor Neurons",
                                                                                                                                               "Cholinergic Interneurons",
                                                                                                                                               "Skeletal Motor Neurons"))

mainlymotor_MFdeep$newnames <- Idents(mainlymotor_MFdeep)
current_ids <- c(0:22)
mainlymotor_MFdeep@meta.data$newnames <- plyr::mapvalues(x = mainlymotor_MFdeep@meta.data$newnames, from = current_ids, to = c("S1",
                                                                                                                               "I1",
                                                                                                                               "S2",
                                                                                                                               "S3",
                                                                                                                               "S3",
                                                                                                                               "V1",
                                                                                                                               "I2",
                                                                                                                               "V2",
                                                                                                                               "I3",
                                                                                                                               "V3",
                                                                                                                               "I4",
                                                                                                                               "V4",
                                                                                                                               "V5",
                                                                                                                               "I5",
                                                                                                                               "V6",
                                                                                                                               "V7",
                                                                                                                               "I6",
                                                                                                                               "V8",
                                                                                                                               "I8",  
                                                                                                                               "V9",
                                                                                                                               "V10",
                                                                                                                               "I7",
                                                                                                                               "S3"))

Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'newnames'

# Recluster data by type of cholinergic neuron
visceral_MNs_MFdeep <- subset(mainlymotor_MFdeep, ident = c(5, 20, 9, 12, 19, 14, 11, 7, 17, 15))
DefaultAssay(visceral_MNs_MFdeep) <- "integrated"
visceral_MNs_MFdeep <- RunPCA(object = visceral_MNs_MFdeep,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(visceral_MNs_MFdeep,ndims=40)
visceral_MNs_MFdeep <- RunUMAP(object = visceral_MNs_MFdeep, reduction = "pca", dims = c(1:15))
visceral_MNs_MFdeep <- FindNeighbors(object = visceral_MNs_MFdeep, reduction = "pca", dims = c(1:15))
visceral_MNs_MFdeep <- FindClusters(visceral_MNs_MFdeep, resolution = 0.6)
visceral_MNs_MFdeep <- RunTSNE(object = visceral_MNs_MFdeep, reduction = "pca", dims = c(1:15))
DefaultAssay(object = visceral_MNs_MFdeep) <- "SCT"
DimPlot(visceral_MNs_MFdeep, split.by = 'expt', pt.size = .2)

skeletal_MNs_MFdeep <- subset(mainlymotor_MFdeep, ident = c(0, 2:4, 22))
DefaultAssay(skeletal_MNs_MFdeep) <- "integrated"
skeletal_MNs_MFdeep <- RunPCA(object = skeletal_MNs_MFdeep,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(skeletal_MNs_MFdeep,ndims=40)
skeletal_MNs_MFdeep <- RunUMAP(object = skeletal_MNs_MFdeep, reduction = "pca", dims = c(1:10))
skeletal_MNs_MFdeep <- FindNeighbors(object = skeletal_MNs_MFdeep, reduction = "pca", dims = c(1:10))
skeletal_MNs_MFdeep <- FindClusters(skeletal_MNs_MFdeep, resolution = 0.6)
skeletal_MNs_MFdeep <- RunTSNE(object = skeletal_MNs_MFdeep, reduction = "pca", dims = c(1:10))
DefaultAssay(object = skeletal_MNs_MFdeep) <- "SCT"
DimPlot(skeletal_MNs_MFdeep, pt.size = .2)

cholIn_MFdeep <- subset(mainlymotor_MFdeep, ident = c(1, 6, 8, 10, 13, 16, 18, 21))
DefaultAssay(cholIn_MFdeep) <- "integrated"
cholIn_MFdeep <- RunPCA(object = cholIn_MFdeep,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(cholIn_MFdeep,ndims=40)
cholIn_MFdeep <- RunUMAP(object = cholIn_MFdeep, reduction = "pca", dims = c(1:15))
cholIn_MFdeep <- FindNeighbors(object = cholIn_MFdeep, reduction = "pca", dims = c(1:15))
cholIn_MFdeep <- FindClusters(cholIn_MFdeep, resolution = 0.6)
cholIn_MFdeep <- RunTSNE(object = cholIn_MFdeep, reduction = "pca", dims = c(1:15))
DefaultAssay(object = cholIn_MFdeep) <- "SCT"
DimPlot(cholIn_MFdeep, split.by = 'expt', pt.size = .2)

#Find markers
visc_markers <- FindAllMarkers(visceral_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
skel_markers <- FindAllMarkers(skeletal_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholin_markers <- FindAllMarkers(cholIn_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
mot_markers <- FindAllMarkers(motor_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'cholinergictypes'
mainmot_markers_bytype <- FindAllMarkers(mainlymotor_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)

# Add skeletal motor neuron subtypes to skeletal_MNs_MFdeep metadata
skeletal_MNs_MFdeep$skeletaltypes <- Idents(skeletal_MNs_MFdeep)
current_ids <- c(0:10)
skeletal_MNs_MFdeep@meta.data$skeletaltypes <- plyr::mapvalues(x = skeletal_MNs_MFdeep@meta.data$skeletaltypes, from = current_ids, to = c("Beta Motor Neurons",
                                                                                                                                           "Gamma Motor Neurons",
                                                                                                                                           "Alpha Motor Neurons",
                                                                                                                                           "Gamma Motor Neurons",
                                                                                                                                           "Gamma Motor Neurons",
                                                                                                                                           "Beta Motor Neurons",
                                                                                                                                           "Alpha Motor Neurons",
                                                                                                                                           "Alpha Motor Neurons",
                                                                                                                                           "Alpha Motor Neurons",
                                                                                                                                           "Gamma Motor Neurons",
                                                                                                                                           "Beta Motor Neurons"))
Idents(object = skeletal_MNs_MFdeep) <- skeletal_MNs_MFdeep@meta.data$'skeletaltypes'
skel_markers_bytype <- FindAllMarkers(skeletal_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)

mainmot_markers_filt <- filter(mainmot_markers_bytype, avg_logFC >= .6, pct.2 <= .3)
mainmot_pct.2_top10 <- mainmot_markers_filt %>% group_by(cluster) %>% top_n(-15, pct.2)
DotPlot(mainlymotor_MFdeep, features = unique(mainmot_pct.2_top10$gene)) + coord_flip() + RotatedAxis()

skel_markers_filt <- filter(skel_markers_bytype, avg_logFC >= .6, pct.2 <= .3)
skel_pct.2_top10 <- skel_markers_filt %>% group_by(cluster) %>% top_n(-10, pct.2)
DotPlot(skeletal_MNs_MFdeep, features = unique(skel_pct.2_top10$gene)) + coord_flip() + RotatedAxis()

# Recluster alpha MNs
alphas_MFdeep <- subset(skeletal_MNs_MFdeep, ident = c(2, 6, 7, 8))
DefaultAssay(alphas_MFdeep) <- "integrated"
alphas_MFdeep <- RunPCA(object = alphas_MFdeep,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(alphas_MFdeep,ndims=40)
alphas_MFdeep <- RunUMAP(object = alphas_MFdeep, reduction = "pca", dims = c(1:15))
alphas_MFdeep <- FindNeighbors(object = alphas_MFdeep, reduction = "pca", dims = c(1:15))
alphas_MFdeep <- FindClusters(alphas_MFdeep, resolution = 0.6)
alphas_MFdeep <- RunTSNE(object = alphas_MFdeep, reduction = "pca", dims = c(1:15))
DefaultAssay(object = alphas_MFdeep) <- "SCT"
DimPlot(alphas_MFdeep, label = F, split.by = 'expt') + NoAxes() 

alpha_markers <- FindAllMarkers(alphas_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
write.csv(alpha_markers %>% group_by(cluster) , file = "alpha_markers_DESeq2.csv")
alpha_markers_filt <- filter(alpha_markers, avg_logFC >= .6, pct.2 <= .3)
alpha_pct.2_top10 <- alpha_markers_filt %>% group_by(cluster) %>% top_n(-20, pct.2)
DotPlot(alphas_MFdeep, features = unique(alpha_pct.2_top10$gene)) + coord_flip() + RotatedAxis()

visc_markers_filt <- filter(visc_markers, avg_logFC >= .6, pct.2 <= .3)
visc_pct.2_top3 <- visc_markers_filt %>% group_by(cluster) %>% top_n(-5, pct.2)
DotPlot(visceral_MNs_MFdeep, features = unique(visc_pct.2_top3$gene), split.by = 'expt', cols = c('red', 'green', 'blue')) + coord_flip() + RotatedAxis()

cholIn_sub_MFdeep <- subset(mainlymotor_MFdeep, ident = c("I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8"))
levels(cholIn_sub_MFdeep) <- c('I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8')
cholin_sub_markers <- FindAllMarkers(cholIn_sub_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholIn_sub_markers_filt <- filter(cholin_sub_markers, avg_logFC >= .6, pct.2 <= .3)
cholIn_sub_pct.2_top3 <- cholIn_sub_markers_filt %>% group_by(cluster) %>% top_n(-3, pct.2)
DotPlot(cholIn_sub_MFdeep, features = unique(cholIn_sub_pct.2_top3$gene)) + coord_flip() + RotatedAxis() # Fig 5a

cholin_markers_filt <- filter(cholin_markers, avg_logFC >= .6, pct.2 <= .3)
cholin_pct.2_top10 <- cholin_markers_filt %>% group_by(cluster) %>% top_n(-10, pct.2)
DotPlot(cholIn_MFdeep, features = unique(cholin_pct.2_top10$gene)) + coord_flip() + RotatedAxis()

# Alpha genes in activity-related GO terms
setwd("/Users/alkaslasimr/Desktop/snRNA seq FACS sorted Feb 2020")
GPCR_activity <- read_excel(file.choose()) # GO term spreadsheet from MGI
GPCRact_symbols <- GPCR_activity$Symbol
ica_in_alphas <- as.data.frame(intersect(ica_symbols, alpha_markers$gene)) %>% mutate(gene = intersect(ica_symbols, alpha_markers$gene))
write.csv(ica_in_alphas , file = "ion channel activity GO genes in alpha markers.csv")
alpha_GO <- alpha_markers[alpha_markers$gene %in% ca_symbols, ]
alpha_GO$GO <- " channel activity"
GABA <- alpha_markers[alpha_markers$gene %in% GABAact_symbols, ]
GABA$GO <- 'GABA receptor activity'
GPCR <- alpha_markers[alpha_markers$gene %in% GPCRact_symbols, ]
GPCR$GO <- 'GPCR  activity'
alpha_GO <- rbind(alpha_GO, GABA)
alpha_GO <- rbind(alpha_GO, GPCR)



## Manuscript figures

# Fig 1
motor_MFdeep$types <- Idents(motor_MFdeep)
current_ids <- c(0:37)
motor_MFdeep@meta.data$types <- plyr::mapvalues(x = motor_MFdeep@meta.data$types, from = current_ids, to = c("Excitatory Neurons", #0
                                                                                                             "Cholinergic Neurons", #1
                                                                                                             "Inhibitory Neurons", #2
                                                                                                             "Inhibitory Neurons", #3
                                                                                                             "Excitatory Neurons", #4
                                                                                                             "Cholinergic Neurons", #5
                                                                                                             "Cholinergic Neurons", #6
                                                                                                             "Cholinergic Neurons", #7
                                                                                                             "Cholinergic Neurons", #8
                                                                                                             "Excitatory Neurons", #9
                                                                                                             "Cholinergic Neurons", #10
                                                                                                             "Excitatory Neurons", #11
                                                                                                             "Inhibitory Neurons", #12
                                                                                                             "Cholinergic Neurons", #13
                                                                                                             "Excitatory Neurons", #14
                                                                                                             "Cholinergic Neurons", #15
                                                                                                             "Excitatory Neurons", #16
                                                                                                             "Cholinergic Neurons", #17
                                                                                                             "Excitatory Neurons", #18
                                                                                                             "Cholinergic Neurons", #19
                                                                                                             "Excitatory Neurons", #20
                                                                                                             "Cholinergic Neurons", #21
                                                                                                             "Cholinergic Neurons", #22
                                                                                                             "Cholinergic Neurons", #23
                                                                                                             "Excitatory Neurons", #24
                                                                                                             "Cholinergic Neurons", #25
                                                                                                             "Excitatory Neurons", #26
                                                                                                             "Excitatory Neurons", #27
                                                                                                             "Excitatory Neurons", #28
                                                                                                             "Excitatory Neurons", #29
                                                                                                             "Excitatory Neurons", #30
                                                                                                             "Cholinergic Neurons", #31
                                                                                                             "Excitatory Neurons", #32
                                                                                                             "Cholinergic Neurons", #33
                                                                                                             "Excitatory Neurons", #34
                                                                                                             "Inhibitory Neurons", #35
                                                                                                             "Excitatory Neurons", #36
                                                                                                             "Excitatory Neurons")) #37
Idents(object = motor_MFdeep) <- motor_MFdeep@meta.data$'types'
levels(motor_MFdeep) <- c('Cholinergic Neurons', 'Excitatory Neurons', 'Inhibitory Neurons')
DimPlot(motor_MFdeep, reduction = 'umap', pt.size = .1, label = T) + NoAxes() # Fig 1a

# Fig 2
mainlymotor_MFdeep$newnames <- Idents(mainlymotor_MFdeep)
current_ids <- c(0:22)
mainlymotor_MFdeep@meta.data$newnames <- plyr::mapvalues(x = mainlymotor_MFdeep@meta.data$newnames, from = current_ids, to = c("S1",
                                                                                                                               "I1",
                                                                                                                               "S2",
                                                                                                                               "S3",
                                                                                                                               "S3",
                                                                                                                               "V1",
                                                                                                                               "I2",
                                                                                                                               "V2",
                                                                                                                               "I3",
                                                                                                                               "V3",
                                                                                                                               "I4",
                                                                                                                               "V4",
                                                                                                                               "V5",
                                                                                                                               "I5",
                                                                                                                               "V6",
                                                                                                                               "V7",
                                                                                                                               "I6",
                                                                                                                               "V8",
                                                                                                                               "I8",  
                                                                                                                               "V9",
                                                                                                                               "V10",
                                                                                                                               "I7",
                                                                                                                               "S3"))

Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'newnames'
colors <- c('lightpink', 'darkgoldenrod3', 'seagreen3', 'cadetblue3', 'lightpink3', 'deeppink1', 'dodgerblue1', 'seagreen2', 'salmon', 'darkorange', 'darkseagreen2', 'darkorchid1', 'indianred1', 'darkturquoise', 'hotpink2', 'darkolivegreen3', 'mediumpurple2', 'maroon', 'plum', 'deeppink2', 'darkslategray2', 'sandybrown')
DimPlot(mainlymotor_MFdeep, reduction = 'umap', label = T, cols = colors) + NoAxes() + NoLegend() # Fig 2a
levels(mainlymotor_MFdeep) <- c('I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'S1', 'S2', 'S3')
VlnPlot(mainlymotor_MFdeep, c('Pax2', 'Zeb2', 'Tns1'), pt.size = 0, ncol = 1, assay= 'SCT', slot = 'scale.data', cols = colors, adjust = 1.8) # Fig 2b
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'cholinergictypes'
mainmot_markers_filt <- filter(mainmot_markers_bytype, avg_logFC >= .6, pct.2 <= .3)
mainmot_pct.2_top10 <- mainmot_markers_filt %>% group_by(cluster) %>% top_n(-15, pct.2)
DotPlot(mainlymotor_MFdeep, features = unique(mainmot_pct.2_top10$gene)) + coord_flip() + RotatedAxis() # Fig 2c
FeaturePlot(mainlymotor_MFdeep, features = c('Slc6a1', 'Fbn2', 'Tns1'), pt.size = .1, ncol = 1) # Fig 2d

# Raw cholinergic neuron dimplot
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'seurat_clusters'
colors <- c('lightpink', 'darkgoldenrod3', 'seagreen3', 'cadetblue3', 'red', 'lightpink3', 'deeppink1', 'dodgerblue1', 'seagreen2', 'salmon', 'darkorange', 'darkseagreen2', 'darkorchid1', 'indianred1', 'darkturquoise', 'hotpink2', 'darkolivegreen3', 'mediumpurple2', 'maroon', 'plum', 'deeppink2', 'darkslategray2', 'sandybrown', 'red')
DimPlot(mainlymotor_MFdeep, reduction = 'umap', label = T, cols = colors) + NoAxes() + NoLegend() # Fig 2a

# Fig 3
skelsub_MFdeep <- subset(mainlymotor_MFdeep, ident = c('S1', 'S2', 'S3'))
colors_skelsub <- c('lightpink', 'seagreen3', 'cadetblue3')
VlnPlot(skelsub_MFdeep, features = c('Rbfox3', 'Esrrg', 'Gfra1'), pt.size = 0, ncol = 1, cols = colors_skelsub) # Fig 3a
aBy_skelsub <- RenameIdents(
  object = skelsub_MFdeep,
  'S1' = 'Beta', 'S2' = 'Alpha', 'S3' = 'Gamma', 'S4' = 'Gamma', 'S5' = 'Gamma')
aBy_skelsub_markers <- FindAllMarkers(aBy_skelsub, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
aBy_skelsub_markers_filt <- filter(aBy_skelsub_markers, avg_logFC >= .6, pct.2 <= .3)
aBy_skelsub_pct.2_top10 <- aBy_skelsub_markers_filt %>% group_by(cluster) %>% top_n(-10, pct.2)
DotPlot(aBy_skelsub, features = unique(aBy_skelsub_pct.2_top10$gene)) + coord_flip() + RotatedAxis() # Fig 3b
VlnPlot(skelsub_MFdeep, features = c('Stk32a', 'Gpr149', 'Nrp2'), pt.size = 0, ncol = 1, cols = colors_skelsub) # Fig 3c

#Fig 4
colors_alpha <- c('darkolivegreen3', 'salmon', 'dodgerblue1', 'plum', 'darkgoldenrod3', 'darkturquoise', 'darkorchid1', 'maroon')
DimPlot(alphas_MFdeep, reduction = 'umap', split.by = 'expt', cols = colors_alpha) + NoAxes()
DimPlot(alphas_MFdeep, reduction = 'umap', label = T, cols = colors_alpha) + NoAxes() + NoLegend()
alpha_markers <- FindAllMarkers(alphas_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
alpha_markers_filt <- filter(alpha_markers, avg_logFC >= .6, pct.2 <= .3)
alpha_pct.2_top10 <- alpha_markers_filt %>% group_by(cluster) %>% top_n(-20, pct.2)
DotPlot(alphas_MFdeep, features = unique(alpha_pct.2_top10$gene)) + coord_flip() + RotatedAxis()

# Fig 5
cholin_markers <- FindAllMarkers(cholIn_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholin_markers_filt <- filter(cholin_markers, avg_logFC >= .6, pct.2 <= .3)
cholin_pct.2_top10 <- cholin_markers_filt %>% group_by(cluster) %>% top_n(-5, pct.2)
DotPlot(cholIn_MFdeep, features = unique(cholin_pct.2_top10$gene)) + coord_flip() + RotatedAxis()
cholIn_sub_MFdeep <- subset(mainlymotor_MFdeep, ident = c("I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8"))
levels(cholIn_sub_MFdeep) <- c('I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8')
cholin_sub_markers <- FindAllMarkers(cholIn_sub_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholIn_sub_markers_filt <- filter(cholin_sub_markers, avg_logFC >= .6, pct.2 <= .3)
cholIn_sub_pct.2_top3 <- cholIn_sub_markers_filt %>% group_by(cluster) %>% top_n(-3, pct.2)
DotPlot(cholIn_sub_MFdeep, features = unique(cholIn_sub_pct.2_top3$gene)) + coord_flip() + RotatedAxis() # Fig 5a
FeaturePlot(cholIn_sub_MFdeep, c('Pitx2', 'Tox'), pt.size = .1, label = T, blend = T, blend.threshold = .1, cols = c('green', 'red')) + NoAxes()
FeaturePlot(cholIn_sub_MFdeep, c('Piezo2', 'Reln'), pt.size = .1)

# Fig 6
colors_visc <- c('lightpink', 'seagreen3', 'mediumpurple2', 'darkgoldenrod3', 'magenta2', 'darkolivegreen3', 'darkturquoise', 'dodgerblue1', 'seagreen2', 'salmon', 'darkorange', 'maroon', 'darkorchid1', 'plum', 'cadetblue3', 'sandybrown')
DimPlot(visceral_MNs_MFdeep, reduction = 'umap', label = T, cols = colors_visc) # Fig 6a
visc_markers <- FindAllMarkers(visceral_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
visc_markers_filt <- filter(visc_markers, avg_logFC >= .6, pct.2 <= .3)
visc_pct.2_top3 <- visc_markers_filt %>% group_by(cluster) %>% top_n(-3, pct.2)
DotPlot(visceral_MNs_MFdeep, features = unique(visc_pct.2_top3$gene)) + coord_flip() + RotatedAxis() # Fig 6b
DimPlot(visceral_MNs_MFdeep, reduction = 'umap', cols = colors_visc, split.by = 'expt', pt.size = .1) # Fig 6c
DotPlot(visceral_MNs_MFdeep, features = c('Gpc3', 'Dach2', 'Sema5a', 'Sst'), split.by = 'expt', cols = c('blue', 'red', 'green')) + coord_flip() + RotatedAxis() # Fig 6b

# Supplementary Fig 3
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'cholinergictypes'
DimPlot(mainlymotor_MFdeep, reduction = 'umap') + NoAxes()
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'newnames'
levels(mainlymotor_MFdeep) <- c('I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'S1', 'S2', 'S3')
VlnPlot(mainlymotor_MFdeep, c('Mpped2', 'Fbn2'), pt.size = 0, ncol = 1, assay= 'SCT', slot = 'scale.data', cols = colors, adjust = 1.8)
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'cholinergictypes'
DimPlot(mainlymotor_MFdeep, reduction = 'umap', split.by = 'expt') + NoAxes() + NoLegend()

# Supplementary Fig 6
FeaturePlot(mainlymotor_MFdeep, 'Prph')

# Supplementary Fig. 7
colors_skeletal <- c('darkturquoise', 'darkgoldenrod3', 'cadetblue3', 'lightpink3', 'dodgerblue1', 'plum', 'seagreen2', 'darkorange', 'darkorchid1', 'indianred1', 'maroon')
DimPlot(skeletal_MNs_MFdeep, label = T, cols = colors_skeletal)
colors_skelsub <- c('lightpink', 'seagreen3', 'cadetblue3')
VlnPlot(skelsub_MFdeep, features = c('Htr1d', 'Atp1a1', 'Wnt7a'), pt.size = .1, ncol = 1, cols = colors_skelsub, adjust = 1.8) 
FeaturePlot(mainlymotor_MFdeep, features = c('Htr1d', 'Atp1a1', 'Wnt7a'), pt.size = .1, ncol = 1) 
FeaturePlot(skeletal_MNs_MFdeep, features = c('Htr1d', 'Atp1a1', 'Wnt7a'), pt.size = .1, ncol = 1) 
VlnPlot(skelsub_MFdeep, c('Sv2b', 'Glis3', 'Rreb1', 'Plekhg1'), pt.size = 0, ncol = 1, assay= 'SCT', slot = 'scale.data', cols = colors_skelsub, adjust = 1.8)

# Supplementary Fig 11
DimPlot(cholIn_MFdeep, reduction = 'umap', label = T)
DimPlot(cholIn_MFdeep, reduction = 'umap', split.by = 'expt') + NoAxes()
cholIn_markers <- FindAllMarkers(cholIn_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholIn_markers_filt <- filter(cholIn_markers, avg_logFC >= .6, pct.2 <= .3)
cholIn_pct.2_top3 <- cholIn_markers_filt %>% group_by(cluster) %>% top_n(-5, pct.2)
DotPlot(cholIn_sub_MFdeep, features = unique(cholIn_pct.2_top3$gene)) + coord_flip() + RotatedAxis() # Fig 5a

# Supplementary Fig 14
FeaturePlot(mainlymotor_MFdeep, features = c('Bnc2', 'Chodl', 'Glis3', 'Gpr149',
                                             'Mpped2', 'Nrp2', 'Piezo2', 'Pitx2',
                                             'Stk32a', 'Sv2b', 'Sv2a', 'Gpc3',
                                             'Dach2', 'Sst', 'Sema5a', 'Reln',
                                             'Tox', 'C1qtnf4', 'Cpne4', 'Grm5',
                                             'Erbb4', 'Kcnq5', 'Plekhg1',
                                             'Rreb1', 'Rbfox3', 'Cacna1e'), pt.size = .1, ncol = 4)

# Supplementary Fig 4
colors <- c('lightpink', 'darkgoldenrod3', 'seagreen3', 'cadetblue3', 'lightpink3', 'deeppink1', 'dodgerblue1', 'seagreen2', 'salmon', 'darkorange', 'darkseagreen2', 'darkorchid1', 'indianred1', 'darkturquoise', 'hotpink2', 'darkolivegreen3', 'mediumpurple2', 'maroon', 'plum', 'deeppink2', 'darkslategray2', 'sandybrown')
DimPlot(mainlymotor_MFdeep, reduction = 'umap', label = T, cols = colors, split.by = 'sex') + NoAxes() + NoLegend()

# Supplementary Table 2
visc_markers <- FindAllMarkers(visceral_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
skel_markers <- FindAllMarkers(skeletal_MNs_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
alpha_markers <- FindAllMarkers(alphas_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
cholin_markers <- FindAllMarkers(cholIn_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
mot_markers <- FindAllMarkers(motor_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'cholinergictypes'
mainmot_markers_bytype <- FindAllMarkers(mainlymotor_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)
Idents(object = mainlymotor_MFdeep) <- mainlymotor_MFdeep@meta.data$'newnames'
mainmot_markers <- FindAllMarkers(mainlymotor_MFdeep, only.pos = TRUE, test.use = "DESeq2", max.cells.per.ident = 100, return.thresh = 0.1, logfc.threshold = 0.5, min.pct = 0.25)



