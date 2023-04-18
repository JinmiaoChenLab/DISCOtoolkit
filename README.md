# DISCOtoolkit 1.0.0

DISCOtoolkit is an R package that allows users to access data and use the tools provided by the [DISCO database](https://www.immunesinglecell.org/). It provides the following functions:

- Filter and download DISCO data based on sample metadata and cell type information
- CELLiD: cell type annotation
- scEnrichment: geneset enrichment using DISCO DEGs
- CellMapper: project data into DISCO atlas

## Requirement

FastIntegration requires the following packages:

-   [R](https://www.r-project.org/) (\>= 4.0.0)
-   [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) (\>= 4.0.0)
-   [SeuratObject](https://cran.r-project.org/web/packages/SeuratObject/index.html) (\>= 4.0.0)
-   [pbmcapply](https://cran.r-project.org/web/packages/pbmcapply/index.html)
-   [stringr](https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html)
-   [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)

## Installation

``` r
devtools::install_github("git@github.com:JinmiaoChenLab/DISCOtoolkit.git")
```

## Basic Usage

### Filter and download DISCO data

``` r
library(DISCOtoolkit)

# find samples from normal lung tissue and sequenced by 10X Genomics platform
# retain samples contain more than 100 Macrophages(or its children)
metadata = FilterDiscoMetadata(
  sample = NULL,
  project = NULL,
  tissue = "lung",
  disease = NULL,
  platform = c("10x3'", "10x5'"),
  sample.type = c("Normal", "Adjacent normal"),
  cell.type = "Macrophage", 
  cell.type.confidence = "high", 
  include.cell.type.children = T, 
  min.cell.per.sample = 50
)
### print information ###
# Retrieving metadata from DISCO database
# Filtering sample
# Retrieving cell type information of each sample from DISCO database
# Retrieving ontology from DISCO database
# 57 samples and 74399 cells were found

# download filtered data into 'disco_data' folder
download.log = DownloadDiscoData(metadata, output.dir = "disco_data")
```

### CELLiD

``` r
library(DISCOtoolkit)
library(Seurat)

metadata = FilterDiscoMetadata(
  sample = "ERX2757110"
)

download.log = DownloadDiscoData(metadata, output.dir = "disco_data")

rna = readRDS("disco_data/ERX2757110.rds")
rna = NormalizeData(rna)
rna.average = AverageExpression(rna)
predict.ct = CELLiDCluster(rna = rna.average$RNA)


# It will download reference data and differential expression gene (DEG) data from DISCO and save them in the 'DISCOtmp' folder by default. You can reuse this data for subsequent CELLiD analyses as follow:

ref.data = readRDS("DISCOtmp/ref_data.rds")
ref.deg = readRDS("DISCOtmp/ref_deg.rds")
predict.ct = CELLiDCluster(rna = rna.average$RNA, ref.data = ref.data, ref.deg = ref.deg)


rna$cell.type = predict.ct$predict_cell_type_1[as.numeric(rna$seurat_clusters)]

DimPlot(rna, group.by = "cell.type", label = T)
```

### scEnrichment

``` r
markers = FindMarkers(rna, ident.1 = 0, only.pos = T, logfc.threshold = 0.5)
cellid.input = data.frame(gene = rownames(markers), logFC = markers$avg_log2FC)
cellid.res = CELLiDEnrichment(cellid.input)

# also it will download 'ref_geneset.rds' in 'DISCOtmp' folder by default,
# You can reuse this data for subsequent CELLiDEnrichment analyses as follow:
ref = readRDS("DISCOtmp/ref_geneset.rds")
cellid.res = CELLiDEnrichment(cellid.input, reference = ref)
```
