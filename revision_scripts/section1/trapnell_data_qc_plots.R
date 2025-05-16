
curr_cell_line <- "A549"

#Read in Trapnell data and convert it to Seurat object
data <- readRDS(paste0(dataDirectory, "raw_data/sciplex_data/GSM4150378_sciPlex3_", curr_cell_line, "_24hrs.RDS"))
rownames(data) <- rowData(data)$gene_short_name
data <- as.Seurat(data, data=NULL)



data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

Idents(data)

VlnPlot(data, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3,pt.size=.1)

plot1 <- FeatureScatter(data, feature1 = "nCount_originalexp", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp")
plot1 + plot2

dim(data)
#Filter cells
nFeature_max <- ifelse(curr_cell_line == "MCF7", 8000, 5000)
data <- subset(data, subset = nFeature_originalexp > 200 & nFeature_originalexp < nFeature_max & percent.mt < 60)

rownames(data)

data@assays$originalexp


data.obj <- CreateSeuratObject(counts=counts, project=paste0(curr_cell_line,"_trapnell"), min.cells = 100)

dim(data)


# Select counts and create new filtered object 
counts <- data@assays$originalexp@counts

data.obj <- CreateSeuratObject(counts=counts, project=paste0(curr_cell_line,"_trapnell"), min.cells = 100)

dim(data.obj)
dim(data)

data.obj[["percent.mt"]] <- PercentageFeatureSet(data.obj, pattern = "^MT-")

Idents(data.obj) <- "a"

VlnPlot(data.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=.1)

plot1 <- FeatureScatter(data.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2







data <- data[rownames(data) %in% rownames(data.obj), colnames(data) %in% colnames(data.obj)]

# data.obj@assays$RNA$counts
data.obj@reductions <- data@reductions
data.obj@meta.data <- data@meta.data
data.obj@reductions$PCA@assay.used <- "RNA"
data.obj@reductions$UMAP@assay.used <- "RNA"
colnames(data.obj@meta.data) <- sapply(colnames(data.obj@meta.data), FUN = function(x) gsub("_originalexp","_RNA", x))

data.obj[["percent.mt"]] <- PercentageFeatureSet(data.obj, pattern = "^MT-")


rownames(data.obj)[grepl("^MT-",rownames(data.obj))]

grepl()

Idents(data.obj) <- data.obj$cell_type

VlnPlot(data.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=.1)

plot1 <- FeatureScatter(data.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Filter cells
nFeature_max <- ifelse(curr_cell_line == "MCF7", 8000, 5000)
data.obj <- subset(data.obj, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_max & percent.mt < 60)

data.obj <- NormalizeData(data.obj)

data.obj <- AddMetaData(data.obj, ifelse(data.obj$dose == 0, "pre", "post"), col.name = "treatment_stage")

data.obj <- FindVariableFeatures(data.obj, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(data.obj)
data.obj <- ScaleData(data.obj, features = all.genes)

data.obj <- RunPCA(data.obj, features = VariableFeatures(object = data.obj))

ElbowPlot(data.obj)

data.obj <- FindNeighbors(data.obj, dims = 1:10)
data.obj <- FindClusters(data.obj, resolution = 0.5)

data.obj <- RunUMAP(data.obj, dims = 1:10)

DimPlot(data.obj, reduction = "umap", group.by = "experiment")