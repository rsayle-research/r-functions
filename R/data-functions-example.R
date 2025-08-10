#' Load and clean data
#' 
#' This function loads a CSV file and performs basic cleaning
#' 
#' @param file_path Path to the CSV file
#' @return A cleaned data frame
#' @examples
#' clean_data("my_data.csv")
clean_data <- function(file_path) {
  data <- read.csv(file_path)
  # Remove rows with all NA values
  data <- data[rowSums(is.na(data)) != ncol(data), ]
  return(data)
}
