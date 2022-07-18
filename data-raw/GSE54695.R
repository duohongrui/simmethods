## Download data
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE54695&format=file&file=GSE54695%5Fdata%5Ftranscript%5Fcounts%2Etxt%2Egz",
              destfile = "./data-raw/GSE54695_data_transcript_counts.txt.gz",
              method = "auto")
## Read data
data <- read.table("./data-raw/GSE54695_data_transcript_counts.txt.gz",
                   header = T,
                   row.names = 1)
## Filter
data <- data[, c(1:80, 161:240)]
## ERCC counts
ERCC_data <- data[grep(rownames(data), pattern = "^ERCC-"), ]
## Delete ERCC counts
data <- data[-grep(rownames(data), pattern = "^ERCC-"), ]
## Sample extra 3950 genes
data <- data[rowSums(data) > 200, ]
set.seed(111)
data <- data[sample(1:nrow(data), 3950), ]
## rbind
data <- as.matrix(round(rbind(data, ERCC_data)))
## group condition (0 for 2i, 1 for serum)
group_condition <- factor(c(rep(0, 80), rep(1, 80)))

## Save
usethis::use_data(data, overwrite = TRUE)
usethis::use_data(group_condition, overwrite = TRUE)
