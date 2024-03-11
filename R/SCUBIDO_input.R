
#' A function which sorts the modern and fossil proxy data sets
#'
#'
#' This function takes the modern data set (modern_data) containing age in BCE/CE in the first
#' column, a modern observational climate timeseries in the second column
#' and then a final array of columns represeting the different xrf elements.
#' The second data set is the fossil data (fossil_data) and this take the same
#' approach as the modern, but contains one less column as there is no
#' climate data included.
#' This function then takes these two data sets and rescales the fossil proxy
#' data to be consistent with the modern xrf data. It then returns a list
#' containing the data required for the next stages of the proxy modelling.
#'
#' @param modern_data a data set containing age in BCE/CE in the first column, a climate data set in the second column and then the third column thereafter contain all of the XRF elements.
#'
#' @param fossil_data a data set containing age in BCE/CE in the first column,the second column and thereafter contain all of the XRF elements used in the modern data
#'
#' @return list of transformed proxy data to be used in the models and a validation dataset
#'
#' @examples
#' # A simple example containing 50 years of modern data and 150 years
#' # of fossil data
#'
#' \dontrun{
#' data(modern_data, fossil_data)
#' SCUBIDO_input(modern_data, fossil_data)
#' }
#' @import dplyr
#' @import magrittr
#' @import utils
#'
#'
#' @export
#'
SCUBIDO_input <- function(modern_data, fossil_data) {

  samplesize <- min(20, nrow(modern_data)) # number of sample to remove
  sample <- head(modern_data, samplesize) # get the data from the modern
  modern_data <- tail(modern_data, -samplesize) # extract this from the modern data

  # Create the validation data frame with the first and last column of the last 20 rows of modern_data
  validation <- data.frame(firstColumn = sample[,1], lastColumn = sample[,ncol(sample)])
  sample <- sample[,-ncol(sample)]   # Remove the last column from the last 20 rows of modern_data

  # get the last rows of modern and put onto fossil
  fossil_data <- rbind(fossil_data, sample)


  # Extract the first column and create age_m dataframe
  time_m <- modern_data %>%
    select(1) %>%
    rename(time_m = 1)

  # Extract the second column and create climate_m dataframe
  temp_m <- modern_data %>%
    select(2) %>%
    rename(temp_m = 1)

  # Extract all columns except the first two and create xrf_m dataframe
  xrf_m <- modern_data %>%
    select(-1, -2)

  # Perform scaling on xrf_m
  xrf_m_resc <- scale(xrf_m)
  xrf_m_means <- attributes(xrf_m_resc)$`scaled:center`
  xrf_m_sds <- attributes(xrf_m_resc)$`scaled:scale`


  # Extract the first column and create age_m dataframe
  time_f <- fossil_data %>%
    select(1) %>%
    rename(time_f = 1)

  # Extract all columns except the first two and create xrf_m dataframe
  xrf_f <- fossil_data %>%
    select(-1,)

  # perform scaling on xrf_f
  xrf_f_resc <- scale(xrf_f, center = xrf_m_means, scale = xrf_m_sds)


  # Create a list to store the results
  sorted_list <- list(
    time_m = as.numeric(as.character(unlist(time_m[[1]]))),
    temp_m = as.numeric(as.character(unlist(temp_m[[1]]))),
    xrf_m_resc = data.matrix(xrf_m_resc),
    time_f = as.numeric(as.character(unlist(time_f[[1]]))),
    xrf_f = data.matrix(xrf_f),
    xrf_f_resc = data.matrix(xrf_f_resc),
    validation = validation
  )

  # Return the result_list (optional)
  class(sorted_list) <- ("sorted")
  return(sorted_list)
}
