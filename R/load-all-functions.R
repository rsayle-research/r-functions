# Load all custom R functions from GitHub
cat("Loading R functions from GitHub...\n")

# Define your function files here - just add new ones to this list
function_files <- c(
  "R/createKASP.R",
  "R/plotting_functions.R",    # Add this when you create it
  "R/stats_functions.R"        # Add this when you create it
  # Just add new files to this list!
)

base_url <- "https://raw.githubusercontent.com/rsayle-research/r-functions/main/"
loaded_functions <- c()

for (file in function_files) {
  tryCatch({
    source(paste0(base_url, file))
    cat("✓ Successfully loaded", basename(file), "\n")
    loaded_functions <- c(loaded_functions, basename(file))
  }, error = function(e) {
    cat("⚠ Skipped", basename(file), "(not found or error)\n")
  })
}

cat("\nLoaded", length(loaded_functions), "function files successfully!\n")
cat("Files loaded:", paste(loaded_functions, collapse = ", "), "\n")
