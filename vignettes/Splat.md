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
    ##             value
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
`prob.group = c(0.3, 0.7)`. In addtion, if we want to simulate three or
more groups, we should obey the rules: \* The length of `prob.group`
vector must always equal to the number of groups. \* The sum of
`prob.group` numeric vector must equal to 1.

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
length of the vector equals to the number of batches. In addtion, if we
want to simulate three or more batches, we should obey the rules: \* The
length of `prob.group` vector always equals to the number of batches. \*
The sum of `prob.group` numeric vector represents the total number of
cells.

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
