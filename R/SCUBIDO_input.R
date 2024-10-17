
#' A function which sorts the modern and fossil proxy data sets
#'
#'
#' This function takes the modern data set (modern_data) containing age in before present (BP) in the first
#' column, a modern observational climate time series in the second column (in anomalies),
#' and then the final columns should contain the data from different xrf elements.
#' The second data set is the fossil data (fossil_data) and this take the same
#' approach as the modern, but contains one less column as there is no
#' climate data included thus age BP in the first and then the remaining columns containing the proxy data.
#' This function then takes these two data sets and re scales the fossil proxy
#' data to be consistent with the modern xrf data. It then returns a list
#' containing the data required for the next stages of the proxy modelling.
#'
#' @param modern_data a data set containing age in BP in the first column, a climate data set in the second column in anomalies, and then the third column thereafter contain all of the XRF elements.
#'
#' @param fossil_data a data set containing age in BP in the first column,the second column and thereafter contain all of the XRF elements used in the modern data
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

# creating our time_grid etc
  time_grid = seq(min(time_m$time_m), max(time_f$time_f), by = 1)
  time_fm <- sort(c(time_f$time_f, time_m$time_m))
  time_all <- sort(unique(c(time_fm, time_grid)))
  duplicate_indices <- which(diff(time_all) == 0)

  # If duplicates are found, issue a warning
  if (length(duplicate_indices) > 0) {
    warning_message <- paste0("Warning: Duplicates found at indices: ",
                              paste(duplicate_indices, collapse = ", please ensure that in your fossil and modern data there are NO overlapping times"))
    stop(warning_message)
  }

  # Create a list to store the results
  sorted_list <- list(
    time_m = as.numeric(as.character(unlist(time_m[[1]]))),
    temp_m = as.numeric(as.character(unlist(temp_m[[1]]))),
    xrf_m_resc = data.matrix(xrf_m_resc),
    time_f = as.numeric(as.character(unlist(time_f[[1]]))),
    xrf_f = data.matrix(xrf_f),
    xrf_f_resc = data.matrix(xrf_f_resc),
    time_grid = time_grid,
    time_fm = time_fm,
    time_all = time_all
  )

  # Return the result_list (optional)
  class(sorted_list) <- ("sorted")
  return(sorted_list)
}
