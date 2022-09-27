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

estimate_result <- simmethods::zinbwaveZinger_estimation(
  ref_data = ref_data,
  other_prior = other_prior,
  verbose = TRUE,
  seed = 111
)

## simulation
other_prior <- list(nCells = 1000,
                    nGenes = 5000,
                    de.prob = 0.2,
                    group.condition = group)
suppressWarnings(
  simulate_result <- simmethods::zinbwaveZinger_simulation(
    ref_data = ref_data,
    parameters = estimate_result[["estimate_result"]],
    other_prior = other_prior,
    return_format = "list",
    verbose = FALSE,
    seed = 111
  )
)

test_that("data size", {
  counts <- simulate_result$simulate_result$count_data
  expect_equal(dim(counts), c(other_prior$nGenes, other_prior$nCells))
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
               other_prior$de.prob)
})
