Splat
================
DUO Hongrui

<!-- github markdown built using 
rmarkdown::render("vignettes/Splat.Rmd", output_format = rmarkdown::github_document())
-->

Here Splat method will be demonstrated clearly and hope that this
document can help you.

## Estimating parameters from a real dataset

Before simulating datasets, it is important to estimate some essential
parameters from a real dataset in order to make the simulated data more
real. If you do not have a single-cell transcriptomics count matrix now,
you can use the data collected in simmethods package by
`simmethods:data` command.

``` r
library(simmethods)
library(SingleCellExperiment)
# Load data
ref_data <- simmethods::data
dim(ref_data)
```

    ## [1] 4000  160

Using `simmethods::Splat_estimation` command to execute the estimation
step.

``` r
estimate_result <- simmethods::Splat_estimation(ref_data = ref_data,
                                                verbose = T,
                                                seed = 10)
```

    ## Estimating parameters using Splat

## Simulating datasets using Splat

After estimating parameter from a real dataset, we will simulate a
dataset based on the learned parameters with different scenarios.

1.  Datasets with default parameters
2.  Determin the number of cells and genes
3.  Simulate two or more groups
4.  Simulate two or more batches
5.  Simulate more groups and batches simutaniously
6.  Return results with different format

### 1. Datasets with default parameters

The reference data contains 160 cells and 4000 genes, if we simulate
datasets with default parameters and then we will obtain a new data
which has the same size as the reference data. In addtion, the simulated
dataset will have one group and one batch of cells.

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "SCE",
                                                seed = 111)
```

    ## nCells: 160 
    ## nGenes: 4000 
    ## nGroups: 1 
    ## de.prob: 0.1 
    ## nBatches: 1

``` r
SCE_result <- simulate_result[["simulate_result"]]
dim(SCE_result)
```

    ## [1] 4000  160

``` r
head(colData(SCE_result))
```

    ## DataFrame with 6 rows and 3 columns
    ##         cell_name       batch       group
    ##       <character> <character> <character>
    ## Cell1       Cell1      Batch1      Group1
    ## Cell2       Cell2      Batch1      Group1
    ## Cell3       Cell3      Batch1      Group1
    ## Cell4       Cell4      Batch1      Group1
    ## Cell5       Cell5      Batch1      Group1
    ## Cell6       Cell6      Batch1      Group1

``` r
head(rowData(SCE_result))
```

    ## DataFrame with 6 rows and 1 column
    ##                 X
    ##       <character>
    ## Gene1       Gene1
    ## Gene2       Gene2
    ## Gene3       Gene3
    ## Gene4       Gene4
    ## Gene5       Gene5
    ## Gene6       Gene6

### 2. Determin the number of cells and genes

In Splat, we can not set `nCells` directly and should set `batchCells`
instead. For example, if we want to simulate 500 cells, we can type
`other_prior = list(batchCells = 500)`. For genes, we can just set
`nGenes`. Here, we simulate a new dataset with 500 cells and 1000 genes:

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "list",
                                                other_prior = list(batchCells = 500,
                                                                   nGenes = 1000),
                                                seed = 111)
```

    ## nCells: 500 
    ## nGenes: 1000 
    ## nGroups: 1 
    ## de.prob: 0.1 
    ## nBatches: 1

``` r
result <- simulate_result[["simulate_result"]][["count_data"]]
dim(result)
```

    ## [1] 1000  500

### 3. Simulate two or more groups

In Splat, we can not set `nGroups` directly and should set `prob.group`
instead. For example, if we want to simulate 2 groups, we can type
`other_prior = list(prob.group = c(0.5, 0.5))`. Note that the sum of
`prob.group` numeric vector must equal to 1, so we can also set
`prob.group = c(0.3, 0.7)`.

In addtion, if we want to simulate three or more groups, we should obey
the rules:

-   The length of `prob.group` vector must always equal to the number of
    groups.
-   The sum of `prob.group` numeric vector must equal to 1.

For demonstration, we will simulate three groups using the learned
parameters.

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "list",
                                                other_prior = list(batchCells = 500,
                                                                   nGenes = 1000,
                                                                   prob.group = c(0.1, 0.3, 0.6)),
                                                seed = 111)
```

    ## nCells: 500 
    ## nGenes: 1000 
    ## nGroups: 3 
    ## de.prob: 0.1 
    ## nBatches: 1

``` r
result <- simulate_result[["simulate_result"]][["count_data"]]
dim(result)
```

    ## [1] 1000  500

``` r
## cell information
cell_info <- simulate_result[["simulate_result"]][["col_meta"]]
table(cell_info$group)
```

    ## 
    ## Group1 Group2 Group3 
    ##     46    156    298

``` r
## gene information
gene_info <- simulate_result[["simulate_result"]][["row_meta"]]
### the proportion of DEGs
table(gene_info$de_gene)[2]/nrow(result) ## de.prob = 0.1
```

    ## yes 
    ## 0.1

We can see that the proportion of differential expressed genes is 0.1
(equals to the default). Next, if we want to know the fold change
between two groups, we will do division with the groups that we are
interested in.

``` r
## fc between group2 and group1
fc_group1_to_group2 <- gene_info$DEFacGroup2/gene_info$DEFacGroup1
## fc between group3 and group1
fc_group1_to_group3 <- gene_info$DEFacGroup3/gene_info$DEFacGroup1
## fc between group3 and group2
fc_group2_to_group3 <- gene_info$DEFacGroup3/gene_info$DEFacGroup2
```

### 4. Simulate two or more batches

In Splat, we can not set `nBatches` directly and should set `batchCells`
instead. For example, if we want to simulate 2 batches, we can type
`other_prior = list(batchCells = c(250, 250))`. Note that the sum of
`batchCells` numeric vector represents the total number of cells and the
length of the vector equals to the number of batches.

In addtion, if we want to simulate three or more batches, we should obey
the rules:

-   The length of `batchCells` vector always equals to the number of
    batches.
-   The sum of `batchCells` numeric vector represents the total number
    of cells.

For demonstration, we will simulate three batches using the learned
parameters.

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "list",
                                                other_prior = list(batchCells = c(200, 300),
                                                                   nGenes = 1000),
                                                seed = 111)
```

    ## nCells: 500 
    ## nGenes: 1000 
    ## nGroups: 1 
    ## de.prob: 0.1 
    ## nBatches: 2

``` r
result <- simulate_result[["simulate_result"]][["count_data"]]
dim(result)
```

    ## [1] 1000  500

``` r
## cell information
cell_info <- simulate_result[["simulate_result"]][["col_meta"]]
table(cell_info$batch)
```

    ## 
    ## Batch1 Batch2 
    ##    200    300

### 5. Simulate more groups and batches simutaniously

As mentioned before, we can set `prob.group` and `batchCells` to
determine the number of groups and batches and we can also set `de.prob`
to specify the proportion of DEGs. Here, we simulate a dataset with
following settings:

-   1000 cells
-   5000 genes
-   three groups with 0.2 proportion of DEGs
-   two batches

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "list",
                                                other_prior = list(batchCells = c(500, 500),
                                                                   nGenes = 5000,
                                                                   de.prob = 0.2,
                                                                   prob.group = c(0.2, 0.3, 0.5)),
                                                seed = 111)
## nCells: 1000 
## nGenes: 5000 
## nGroups: 3 
## de.prob: 0.2 
## nBatches: 2
result <- simulate_result[["simulate_result"]][["count_data"]]
dim(result)
## [1] 5000 1000
## cell information
cell_info <- simulate_result[["simulate_result"]][["col_meta"]]
table(cell_info$batch)
## 
## Batch1 Batch2 
##    500    500
table(cell_info$group)
## 
## Group1 Group2 Group3 
##    186    321    493
## gene information
gene_info <- simulate_result[["simulate_result"]][["row_meta"]]
### proportion of DEGs
table(gene_info$de_gene)[2]/nrow(result)
##    yes 
## 0.1888
### fc
fc_group2_to_group3 <- gene_info$DEFacGroup3/gene_info$DEFacGroup2
```

### 6. Return results with different format

In simmethods package, we provide four formats of results to users
without data format conversion, including `list`,
`SingleCellExperiment`, `Seurat` and `h5ad`. The previous three formats
are compatible with R environment and the last `h5ad` format is suitable
for **Python** environment and can be imported by `scanpy.read_h5ad`
function.

#### list

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "list",
                                                other_prior = list(batchCells = c(100, 100),
                                                                   nGenes = 1000,
                                                                   de.prob = 0.1,
                                                                   prob.group = c(0.2, 0.3, 0.5)),
                                                seed = 111)
## nCells: 200 
## nGenes: 1000 
## nGroups: 3 
## de.prob: 0.1 
## nBatches: 2
str(simulate_result)
## List of 2
##  $ simulate_result   :List of 3
##   ..$ count_data: int [1:1000, 1:200] 14 19 21 42 0 16 17 708 20 17 ...
##   .. ..- attr(*, "dimnames")=List of 2
##   .. .. ..$ : chr [1:1000] "Gene1" "Gene2" "Gene3" "Gene4" ...
##   .. .. ..$ : chr [1:200] "Cell1" "Cell2" "Cell3" "Cell4" ...
##   ..$ col_meta  :'data.frame':   200 obs. of  3 variables:
##   .. ..$ cell_name: chr [1:200] "Cell1" "Cell2" "Cell3" "Cell4" ...
##   .. ..$ batch    : chr [1:200] "Batch1" "Batch1" "Batch1" "Batch1" ...
##   .. ..$ group    : Factor w/ 3 levels "Group1","Group2",..: 2 2 3 2 3 3 3 2 3 3 ...
##   ..$ row_meta  :'data.frame':   1000 obs. of  7 variables:
##   .. ..$ gene_name     : chr [1:1000] "Gene1" "Gene2" "Gene3" "Gene4" ...
##   .. ..$ de_gene       : chr [1:1000] "no" "no" "no" "no" ...
##   .. ..$ BatchFacBatch1: num [1:1000] 0.91 0.977 1.054 1.171 1.002 ...
##   .. ..$ BatchFacBatch2: num [1:1000] 1.105 0.95 0.76 0.776 0.996 ...
##   .. ..$ DEFacGroup1   : num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
##   .. ..$ DEFacGroup2   : num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
##   .. ..$ DEFacGroup3   : num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
##  $ simulate_detection:'data.frame':  1 obs. of  4 variables:
##   ..$ Function_Call     : chr "simulate_result<-splatter::splatSimulate(parameters,method=submethod,verbose=verbose)"
##   ..$ Elapsed_Time_sec  : num 0.55
##   ..$ Total_RAM_Used_MiB: num 6.9
##   ..$ Peak_RAM_Used_MiB : num 36.8
counts <- simulate_result[["simulate_result"]][["count_data"]]
## cell information
cell_info <- simulate_result[["simulate_result"]][["col_meta"]]
head(cell_info)
##       cell_name  batch  group
## Cell1     Cell1 Batch1 Group2
## Cell2     Cell2 Batch1 Group2
## Cell3     Cell3 Batch1 Group3
## Cell4     Cell4 Batch1 Group2
## Cell5     Cell5 Batch1 Group3
## Cell6     Cell6 Batch1 Group3
## gene information
gene_info <- simulate_result[["simulate_result"]][["row_meta"]]
head(gene_info)
##       gene_name de_gene BatchFacBatch1 BatchFacBatch2 DEFacGroup1 DEFacGroup2 DEFacGroup3
## Gene1     Gene1      no      0.9098860      1.1054169           1           1           1
## Gene2     Gene2      no      0.9774161      0.9501320           1           1           1
## Gene3     Gene3      no      1.0541276      0.7597880           1           1           1
## Gene4     Gene4      no      1.1708139      0.7762219           1           1           1
## Gene5     Gene5      no      1.0017116      0.9963063           1           1           1
## Gene6     Gene6      no      0.7654326      1.1497335           1           1           1
```

#### SingleCellExperiment

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "SingleCellExperiment",
                                                other_prior = list(batchCells = c(100, 100),
                                                                   nGenes = 1000,
                                                                   de.prob = 0.1,
                                                                   prob.group = c(0.2, 0.3, 0.5)),
                                                seed = 111)
## nCells: 200 
## nGenes: 1000 
## nGroups: 3 
## de.prob: 0.1 
## nBatches: 2
counts <- counts(simulate_result[["simulate_result"]])
## cell information
cell_info <- as.data.frame(colData(simulate_result[["simulate_result"]]))
head(cell_info)
##       cell_name  batch  group
## Cell1     Cell1 Batch1 Group2
## Cell2     Cell2 Batch1 Group2
## Cell3     Cell3 Batch1 Group3
## Cell4     Cell4 Batch1 Group2
## Cell5     Cell5 Batch1 Group3
## Cell6     Cell6 Batch1 Group3
## gene information
gene_info <- as.data.frame(rowData(simulate_result[["simulate_result"]]))
head(gene_info)
##       gene_name de_gene BatchFacBatch1 BatchFacBatch2 DEFacGroup1 DEFacGroup2 DEFacGroup3
## Gene1     Gene1      no      0.9098860      1.1054169           1           1           1
## Gene2     Gene2      no      0.9774161      0.9501320           1           1           1
## Gene3     Gene3      no      1.0541276      0.7597880           1           1           1
## Gene4     Gene4      no      1.1708139      0.7762219           1           1           1
## Gene5     Gene5      no      1.0017116      0.9963063           1           1           1
## Gene6     Gene6      no      0.7654326      1.1497335           1           1           1
```

#### Seurat

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "Seurat",
                                                other_prior = list(batchCells = c(100, 100),
                                                                   nGenes = 1000,
                                                                   de.prob = 0.1,
                                                                   prob.group = c(0.2, 0.3, 0.5)),
                                                seed = 111)
## nCells: 200 
## nGenes: 1000 
## nGroups: 3 
## de.prob: 0.1 
## nBatches: 2
seurat_result <- simulate_result[["simulate_result"]]
## Overview
seurat_result
## An object of class Seurat 
## 1000 features across 200 samples within 1 assay 
## Active assay: originalexp (1000 features, 0 variable features)
## count matrix
counts <- seurat_result@assays$originalexp@counts
counts[1:10, 1:10]
## 10 x 10 sparse Matrix of class "dgCMatrix"
##    [[ suppressing 10 column names 'Cell1', 'Cell2', 'Cell3' ... ]]
##                                            
## Gene1   14  .  10   1  1  .  .  25  261   1
## Gene2   19  .  18   7  .  1  .  22  265   .
## Gene3   21  2   6   5  .  .  1   6  136   1
## Gene4   42  .   8  16  .  1  4  41  603   2
## Gene5    .  .   3   3  1  .  .   8  127   .
## Gene6   16  .   7  10  3  .  .  29  421   1
## Gene7   17  .   8   5  .  .  .  21  110   6
## Gene8  708 13 565 265 52 37 16 686 7368 101
## Gene9   20  .  24   4  7  .  .  27  199   1
## Gene10  17  1  17   7  .  .  .  20  206   4
## cell information
cell_info <- seurat_result@meta.data
head(cell_info)
##          orig.ident nCount_originalexp nFeature_originalexp cell_name  batch  group
## Cell1 SeuratProject              65961                  977     Cell1 Batch1 Group2
## Cell2 SeuratProject               1140                  296     Cell2 Batch1 Group2
## Cell3 SeuratProject              45613                  954     Cell3 Batch1 Group3
## Cell4 SeuratProject              26289                  905     Cell4 Batch1 Group2
## Cell5 SeuratProject               3266                  497     Cell5 Batch1 Group3
## Cell6 SeuratProject               2002                  389     Cell6 Batch1 Group3
## gene information
gene_info <- seurat_result@assays[["originalexp"]]@meta.features
head(gene_info)
##       gene_name de_gene BatchFacBatch1 BatchFacBatch2 DEFacGroup1 DEFacGroup2 DEFacGroup3
## Gene1     Gene1      no      0.9098860      1.1054169           1           1           1
## Gene2     Gene2      no      0.9774161      0.9501320           1           1           1
## Gene3     Gene3      no      1.0541276      0.7597880           1           1           1
## Gene4     Gene4      no      1.1708139      0.7762219           1           1           1
## Gene5     Gene5      no      1.0017116      0.9963063           1           1           1
## Gene6     Gene6      no      0.7654326      1.1497335           1           1           1
```

#### h5ad

If we select `h5ad` format, it is not possible to return the result in
R, so you can get the path where the `h5ad` files save to and we can go
to the path and read it in **Python** by `scanpy.read_h5ad` function (if
you have already installed Python and **scanpy** module).

``` r
simulate_result <- simmethods::Splat_simulation(parameters = estimate_result[["estimate_result"]],
                                                return_format = "h5ad",
                                                other_prior = list(batchCells = c(100, 100),
                                                                   nGenes = 1000,
                                                                   de.prob = 0.1,
                                                                   prob.group = c(0.2, 0.3, 0.5)),
                                                seed = 111)
## nCells: 200 
## nGenes: 1000 
## nGroups: 3 
## de.prob: 0.1 
## nBatches: 2
## Creating h5Seurat file for version 3.1.5.9900
## Adding counts for originalexp
## Adding data for originalexp
## No variable features found for originalexp
## Adding feature-level metadata for originalexp
## Validating h5Seurat file
## Adding data from originalexp as X
## Transfering meta.features to var
## Adding counts from originalexp as raw
## Transfering meta.features to raw/var
## Transfering meta.data to obs
## Your data has been save to C:/Users/duoho/AppData/Local/Temp/RtmpGeYdMc/20220716133958.h5ad
save_path <- simulate_result[["simulate_result"]][["save_path"]]
save_path
## [1] "C:/Users/duoho/AppData/Local/Temp/RtmpGeYdMc/20220716133958.h5ad"
```

Now, we can go to the path and check the data. Here, we read the `h5ad`
file in R using **reticulate** R package (note that **Python** and
**scanpy** module must have been installed).

``` r
## install.packages("reticulate")
scanpy <- reticulate::import("scanpy")
data <- scanpy$read_h5ad(save_path)
data ## Read h5ad file successfully
## AnnData object with n_obs × n_vars = 200 × 1000
##     obs: 'orig.ident', 'nCount_originalexp', 'nFeature_originalexp', 'cell_name', 'batch', 'group'
##     var: 'gene_name', 'de_gene', 'BatchFacBatch1', 'BatchFacBatch2', 'DEFacGroup1', 'DEFacGroup2', 'DEFacGroup3'
```
