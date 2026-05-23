reactome <- read.table(
  "data/Ensembl2Reactome_All_Levels.human.txt",
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

geneiddata <- read.table(
  "data/gProfiler_hsapiens_2026-05-06_12-58-07.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

colnames(reactome)[colnames(reactome) == "V1"] <- "gene"
colnames(reactome)[colnames(reactome) == "V2"] <- "term"
colnames(reactome)[colnames(reactome) == "V3"] <- "description"

