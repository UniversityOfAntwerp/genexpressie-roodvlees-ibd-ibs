# List of required packages
required_packages <- c(
  "GEOquery",   # GEO data import
  "Biobase",    # ExpressionSet handling: exprs(), pData()
  "lme4",       # Linear mixed models
  "lmerTest",   # p-values for lmer() models
  "readr",      # Reading/writing files
)

# Check if packages are installed
missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

# Stop with clear message if anything is missing
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following packages are missing:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall CRAN packages with:\n",
      'install.packages(c("lme4", "lmerTest", "ggplot2", "dplyr", "tibble", "readr", "stringr"))',
      "\n\nInstall Bioconductor packages with:\n",
      'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")',
      "\n",
      'BiocManager::install(c("GEOquery", "Biobase", "limma"))'
    )
  )
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

message("All required libraries loaded successfully.")
