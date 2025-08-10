# Load all custom R functions
# Source this file to load all your custom functions at once

# Source all R files in the R/ directory
r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  source(file)
}

cat("All custom R functions loaded successfully!\n")
cat("Available functions:\n")
cat("- clean_data() - Basic data cleaning\n")
cat("- createKASP() - Generate KASP sequences from variant data\n")
# Add more function names as you create them
