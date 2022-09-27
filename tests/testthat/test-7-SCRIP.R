data <- scater::mockSCE()
ref_data <- as.matrix(SingleCellExperiment::counts(data))

## estimation
estimate_result <- simmethods::SCRIP_estimation(
  ref_data = ref_data,
  verbose = FALSE,
  seed = 111
)

## simulation
mode <- c("GP-trendedBCV",
          "GP-commonBCV",
          "BGP-commonBCV",
          "BP",
          "BGP-trendedBCV")

for(mode in mode){
  other_prior <- list(prob.group = c(0.4, 0.3, 0.3),
                      batchCells = c(80, 80),
                      de.prob = 0.2,
                      nGenes = 500,
                      mode = mode)
  suppressWarnings(
    simulate_result <- simmethods::SCRIP_simulation(
      ref_data = ref_data,
      parameters = estimate_result[["estimate_result"]],
      other_prior = other_prior,
      return_format = "list",
      verbose = FALSE,
      seed = 111
    )
  )
  print(mode)
}

test_that("the class of estimation result", {
  expect_equal(as.character(class(estimate_result$estimate_result)),
               "SplatParams")
})

test_that("cell information", {
  cell_info <- simulate_result$simulate_result$col_meta
  counts <- simulate_result$simulate_result$count_data
  expect_equal(colnames(cell_info), c("cell_name", "batch", "group"))
  expect_equal(dim(counts), c(other_prior$nGenes, sum(other_prior$batchCells)))
})

test_that("gene information", {
  gene_info <- simulate_result$simulate_result$row_meta
  ## proportion of DEGs
  expect_equal(round(sum(gene_info$de_gene == "yes")/nrow(gene_info), digits = 1),
               other_prior$de.prob)
})

