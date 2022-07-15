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
