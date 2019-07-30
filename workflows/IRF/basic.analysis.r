# load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(doParallel)
registerDoParallel(cores=8)

# load data
data1<-Read10X(data.dir = "1/")
pbmc1 <- CreateSeuratObject(counts = data1, project = "IRF1", min.cells = 3, min.features = 200)
data2<-Read10X(data.dir = "2/")
pbmc2 <- CreateSeuratObject(counts = data2, project = "IRF2", min.cells = 3, min.features = 200)
data3<-Read10X(data.dir = "3/")
pbmc3 <- CreateSeuratObject(counts = data3, project = "IRF3", min.cells = 3, min.features = 200)
data4<-Read10X(data.dir = "4/")
pbmc4 <- CreateSeuratObject(counts = data4, project = "IRF4", min.cells = 3, min.features = 200)
data5<-Read10X(data.dir = "5/")
pbmc5 <- CreateSeuratObject(counts = data5, project = "IRF5", min.cells = 3, min.features = 200)
data6<-Read10X(data.dir = "6/")
pbmc6 <- CreateSeuratObject(counts = data6, project = "IRF6", min.cells = 3, min.features = 200)
data7<-Read10X(data.dir = "7/")
pbmc7 <- CreateSeuratObject(counts = data7, project = "IRF7", min.cells = 3, min.features = 200)
data8<-Read10X(data.dir = "8/")
pbmc8 <- CreateSeuratObject(counts = data8, project = "IRF8", min.cells = 3, min.features = 200)
data9<-Read10X(data.dir = "9/")
pbmc9 <- CreateSeuratObject(counts = data9, project = "IRF9", min.cells = 3, min.features = 200)
data10<-Read10X(data.dir = "10/")
pbmc10 <- CreateSeuratObject(counts = data10, project = "IRF10", min.cells = 3, min.features = 200)
data11<-Read10X(data.dir = "11/")
pbmc11 <- CreateSeuratObject(counts = data11, project = "IRF11", min.cells = 3, min.features = 200)
data12<-Read10X(data.dir = "12/")
pbmc12 <- CreateSeuratObject(counts = data12, project = "IRF12", min.cells = 3, min.features = 200)

# merge data
pbmc<- merge(pbmc1, y = c(pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10,pbmc11,pbmc12), add.cell.ids = c("1","2","3","4","5","6","7","8","9","10","11","12"), project = "IRF")

# save raw data
save(pbmc,file="raw.Robj")

# check cell number in samples
table(pbmc@meta.data$orig.ident)


# calculate the mito gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# relevel
levels(pbmc)<-c("IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9","IRF10","IRF11","IRF12")

# before qc plots-violin
pdf("plot.1.qc.before.pdf",paper = "a4r", width = 0, height = 0)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0)
dev.off()

# before qc plots-cor
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("plot.1.qc.before.cor.1.pdf",paper = "a4r", width = 0, height = 0)
plot1
dev.off()
pdf("plot.1.qc.before.cor.2.pdf",paper = "a4r", width = 0, height = 0)
plot2
dev.off()

# qc
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 50)

# check cell numbers after qc
table(pbmc@meta.data$orig.ident)

# after qc plots-violin
pdf("plot.1.qc.after.pdf",paper = "a4r", width = 0, height = 0)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0)
dev.off()

# after qc plots-cor
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("plot.1.qc.after.cor.1.pdf",paper = "a4r", width = 0, height = 0)
plot1
dev.off()
pdf("plot.1.qc.after.cor.2.pdf",paper = "a4r", width = 0, height = 0)
plot2
dev.off()

# normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# find variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10

# variable genes plots
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("plot.2.vgenes.1.pdf",paper = "a4r", width = 0, height = 0)
plot1
dev.off()
pdf("plot.2.vgenes.2.pdf",paper = "a4r", width = 0, height = 0)
plot2
dev.off()

# scale
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Identifying the true dimensionality of a dataset 
pbmc <- JackStraw(pbmc, num.replicate = 100)

# cluster
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.3)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- FindClusters(pbmc, resolution = 0.1)
Idents(pbmc)<-pbmc@meta.data$RNA_snn_res.1
table(Idents(pbmc))

# umap
pbmc <- RunUMAP(pbmc, dims = 1:30)
pdf("plot.4.umap.clusters.pdf",paper = "a4r", width = 0, height = 0)
DimPlot(pbmc, reduction = "umap")
dev.off()

# save obj
saveRDS(pbmc, file = "clusters.rds")

# find marker genes
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

