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

suppressWarnings(
  estimate_result <- simmethods::scDesign2_estimation(
    ref_data = ref_data,
    other_prior = other_prior,
    verbose = TRUE,
    seed = 111
  )
)

## simulation
test_that("group number error", {
  expect_error({
    simulate_result <- simmethods::scDesign2_simulation(
      parameters = estimate_result[["estimate_result"]],
      other_prior = list(nCells = 500,
                         prob.group = c(0.2, 0.3, 0.5)),
      return_format = "list",
      verbose = FALSE,
      seed = 111
    )
  }, info = "Cell type proportion should have the same length as the number of models.")
})

simulate_result <- simmethods::scDesign2_simulation(
  parameters = estimate_result[["estimate_result"]],
  other_prior = list(nCells = 500,
                     prob.group = c(0.5, 0.5)),
  return_format = "list",
  verbose = FALSE,
  seed = 111
)

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(ncol(counts), 500)
})

test_that("cell information", {
  cell_info <- simulate_result$simulate_result$col_meta
  ## proportion of DEGs
  expect_equal(as.numeric(table(cell_info$group)/nrow(cell_info)),
               c(0.5, 0.5))
})

# test_that("gene information", {
#   gene_info <- simulate_result$simulate_result$row_meta
#   ## proportion of DEGs
#   expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
#                other_prior$de.prob)
# })
