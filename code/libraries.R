# List of required packages
required_packages <- c(
  "GEOquery",   # GEO data import
  "Biobase",    # ExpressionSet handling: exprs(), pData()
  "lme4",       # Linear mixed models
  "lmerTest",   # p-values for lmer() models
  "readr",      # Reading/writing files
)

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

message("All required libraries loaded successfully.")
