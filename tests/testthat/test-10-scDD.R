data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))
group <- as.data.frame(SingleCellExperiment::colData(data))
group <- as.numeric(as.factor(group$Treatment))

test_that("group information", {
  expect_equal(class(group), "numeric")
  expect_equal(min(group), 1)
})

## estimation
other_prior <- list(group.condition = group)

test_that("group information exists", {
  expect_true(!is.null(other_prior$group.condition))
})

estimate_result <- simmethods::scDD_estimation(
  ref_data = ref_data,
  other_prior = other_prior,
  verbose = TRUE,
  seed = 111
)

## simulation
other_prior <- list(nCells = 1000)
suppressWarnings(
  simulate_result <- simmethods::scDD_simulation(
    parameters = estimate_result[["estimate_result"]],
    other_prior = other_prior,
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("the class of estimation result", {
  expect_equal(as.character(class(estimate_result$estimate_result)),
               "SCDDParams")
})

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(ncol(counts), other_prior$nCells)
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(unique(gene_info$DEstatus),
               c("EP", "EE"))
})
