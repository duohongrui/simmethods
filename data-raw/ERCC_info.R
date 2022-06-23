download.file(url = "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt",
              destfile = "./data-raw/ERCC_info.txt",
              method = "auto")
ERCC_info <- read.table("./data-raw/ERCC_info.txt", header = TRUE, row.names = 1, sep = "\t")
colnames(ERCC_info) <- c("ERCC_id",
                         "subgroup",
                         "con_Mix1_attomoles_ul",
                         "con_Mix2_attomoles_ul",
                         "expected_fc",
                         "log2_fc")
usethis::use_data(ERCC_info, overwrite = TRUE)
