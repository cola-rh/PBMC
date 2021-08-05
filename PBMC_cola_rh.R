setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/PBMC")
library(cola)
library(scRNAseq)

library(Seurat)
pbmc.data = Read10X(data.dir = "filtered_gene_bc_matrices")
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/PBMC/PBMC_data.rds')

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(pbmc)
data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/PBMC/PBMC_data.rds')
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

Seurat_class = as.character(Idents(pbmc))

data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/PBMC/PBMC_data.rds')
mat = adjust_matrix(mat)

rh = hierarchical_partition(mat, subset = 500, cores = 4, 
    anno = data.frame(Seurat_class = Seurat_class))
saveRDS(rh, file = "PBMC_cola_rh.rds")


cola_report(rh, output = "PBMC_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'PBMC'")
