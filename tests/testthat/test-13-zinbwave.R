data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))

## estimation
suppressWarnings(
  estimate_result <- simmethods::zinbwave_estimation(
    ref_data = ref_data,
    verbose = TRUE,
    seed = 111
  )
)

## simulation
suppressWarnings(
  simulate_result <- simmethods::zinbwave_simulation(
    parameters = estimate_result[["estimate_result"]],
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("the class of estimation result", {
  expect_equal(as.character(class(estimate_result$estimate_result)),
               "ZINBParams")
})

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(dim(counts), dim(ref_data))
})

