ref_data <- simmethods::data

test_that("ERCC", {
  expect_true(!is.null(grep(rownames(ref_data), pattern = "^ERCC-")))
})

## estimation
suppressWarnings(
  estimate_result <- powsimR_estimation(
    ref_data = ref_data,
    other_prior = list(RNAseq = "singlecell",
                       Protocol = "UMI",
                       Normalisation = "scran",
                       dilution.factor = 50000,
                       volume = 1),
    verbose = TRUE,
    seed = 111
  )
)

## simulation
other_prior <- list(prob.group = c(0.4, 0.6),
                    nCells = 1000,
                    de.prob = 0.2,
                    nGenes = 500)
suppressWarnings(
  simulate_result <- simmethods::powsimR_simulation(
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
  expect_equal(dim(counts), c(other_prior$nGenes, sum(other_prior$nCells)))
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
               other_prior$de.prob)
})

