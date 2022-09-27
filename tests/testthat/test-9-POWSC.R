data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))

## estimation
estimate_result <- simmethods::POWSC_estimation(
  ref_data = ref_data,
  verbose = FALSE,
  seed = 111
)

## simulation
other_prior <- list(nCells = 1000,
                    de.prob = 0.2)
suppressWarnings(
  simulate_result <- simmethods::POWSC_simulation(
    parameters = estimate_result[["estimate_result"]],
    other_prior = other_prior,
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("cell information", {
  cell_info <- simulate_result$simulate_result$col_meta
  counts <- simulate_result$simulate_result$count_data
  expect_equal(colnames(cell_info), c("cell_name", "group"))
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
               other_prior$de.prob)
})

