data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))

## estimation
estimate_result <- simmethods::Simple_estimation(
  ref_data = ref_data,
  verbose = FALSE,
  seed = 111
)

## simulation
other_prior <- list(nCells = 500,
                    nGenes = 1000)
suppressWarnings(
  simulate_result <- simmethods::Simple_simulation(
    parameters = estimate_result[["estimate_result"]],
    other_prior = other_prior,
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("the class of estimation result", {
  expect_equal(as.character(class(estimate_result$estimate_result)),
               "SimpleParams")
})

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(dim(counts), c(other_prior$nGenes, other_prior$nCells))
})

