data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))

## simulation
other_prior <- list(nCells = 500,
                    nGroups = 3,
                    prob.group = c(0.2, 0.3, 0.5),
                    de.prob = 0.3)
suppressWarnings(
  simulate_result <- simmethods::scDesign_simulation(
    ref_data = ref_data,
    other_prior = other_prior,
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(ncol(counts), other_prior$nCells)
})

test_that("cell information", {
  cell_info <- simulate_result$simulate_result$col_meta
  ## proportion of DEGs
  expect_equal(as.numeric(table(cell_info$group)/nrow(cell_info)),
               other_prior$prob.group)
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
               other_prior$de.prob)
})
