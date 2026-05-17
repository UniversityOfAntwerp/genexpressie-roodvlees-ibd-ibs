required_packages <- c(
  "GEOquery",   
  "Biobase",       
  "lmerTest",  
  "readr"      
)

invisible(lapply(required_packages, library, character.only = TRUE))

