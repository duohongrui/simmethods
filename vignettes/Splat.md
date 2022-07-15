Splat
================
DUO Hongrui

<!-- github markdown built using 
rmarkdown::render("vignettes/getting_started.Rmd", output_format = rmarkdown::github_document())
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

### Simulate two or more groups
