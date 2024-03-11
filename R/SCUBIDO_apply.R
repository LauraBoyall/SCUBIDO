#' Apply the relationship between modern climate and XRF core scanning data and apply to the fossil data
#'
#'
#' The SCUBIDO_apply function takes what is learnt in the modern observational
#'     period in \code{\link{SCUBIDO_cal}} and applies this relationship to the fossil data.
#'
#'
#' @param sorted the sorted data set coming from the first function
#' @param calibration_data is the output list from the previous SCUBIDO_cal function.
#' @param temp_grid a sequence range between two values for the model to pick the temperature from. if working with anomalies use seq(-3, 3, length = 50)
#' @param plot_graph a function to plot the output of the SCUBIDO_apply function
#' @param print a function which prints out the output
#'
#' @import ggplot2
#' @import R2jags
#' @import mvtnorm
#' @import fitdistrplus
#' @import mvtnorm
#' @import ggpubr
#'
#'
#' @return a data set containing the output
#'
#'
#' @examples
#' \dontrun{
#' temp_grid <- seq(3,3, length = 50)
#' MDP <- SCUBIDO_apply(calibration_data, temp_grid, plot_graph = T, print = F)
#'}
#'
#'
#SCUBIDO_apply <- function(calibration,sorted, temp_grid, plot_graph = FALSE, print = FALSE){ #technically I don't think we need the sorted here as it should be stored in calibration
SCUBIDO_apply <- function(calibration, temp_grid, plot_graph = FALSE, print = FALSE){
    UseMethod("SCUBIDO_apply")
}

#' @export
SCUBIDO_apply.calibration <- function(calibration,temp_grid, plot_graph = FALSE, print = FALSE) {
#SCUBIDO_apply.calibration.sorted <- function(calibration,sorted,temp_grid, plot_graph = FALSE, print = FALSE) { #technically I don't think we need the sorted here as it should be stored in calibration
  MDP_file <- paste0("MDPs_mising", Sys.Date() , ".rds")
  message("Warning: This operation may take a significant amount of time to complete.
  It is working in the background and you can see printed in the console counting down through
  the number of data points. We advise leaving this code running in background jobs until the counting
  down has stopped.")
  cal_run_results <- calibration$sims.list
  N_rows_f <- nrow(calibration$sorted$xrf_f_resc)
  temp_grid <- seq(-3, 3, length = 50)

  if (!file.exists(MDP_file)) {
    MDP_raw <- MDP <- matrix(0, nrow = N_rows_f, ncol = length(temp_grid))


    for (i in 1:N_rows_f) {
      print(N_rows_f - i)
      for (k in 1:length(temp_grid)) {
        for (j in 1:100) {
          beta0 <- cal_run_results$beta0[j, ]
          beta1 <- cal_run_results$beta1[j, ]
          beta2 <- cal_run_results$beta2[j, ]
          Sigma_inv <- cal_run_results$Sigma_inv[j, , ]
          Sigma <- solve(Sigma_inv)
          mean <- beta0 + beta1 * temp_grid[k] + beta2 * (temp_grid[k]^2)
          MDP_raw[i, k] <- MDP_raw[i, k] +
            dmvnorm(calibration$sorted$xrf_f_resc[i, ],
                    mean = mean,
                    sigma = Sigma,
                    log = TRUE
            )
        }
      }
      c <- max(MDP_raw[i, ])
      MDP[i, ] <- exp((MDP_raw[i, ] - c) / nrow(cal_run_results$beta0))
    }
    MDP_final <- sweep(MDP, 1, apply(MDP, 1, "sum"), "/")
    saveRDS(MDP_final, file = MDP_file)
  } else {
    MDP_final <- readRDS(file = MDP_file)
  }

  MDP_mean <- round(MDP_final, 5) %*% temp_grid
  MDP_sd <- sqrt(abs(round(MDP_final, 5) %*% (temp_grid^2) - MDP_mean^2))
  MDP_prec <- 1 / (MDP_sd^2)
  MDP_high <- MDP_mean + 2 * MDP_sd
  MDP_low <- MDP_mean - 2 * MDP_sd
  MDP_prec[is.infinite(MDP_prec)] <- 1e10

  MDPs <- list(
    time_f = calibration$sorted$time_f,
    MDP_mean = MDP_mean,
    MDP_low = MDP_low,
    MDP_high = MDP_high,
    MDP_prec = MDP_prec,
    sd_rw = cal_run_results,
    sorted = calibration$sorted,
    calibration = calibration$sims.list
  )

  # Assign the calibration$sorted to a variable in the global environment
  #assign("MDP_list", MDP_list, envir = .GlobalEnv)

  df <- data.frame(
    time_f = calibration$sorted$time_f,  # Assuming 'time_f' is the correct column name
    MDP_mean = MDP_mean,
    MDP_low = MDP_low,
    MDP_high = MDP_high,
    MDP_prec = MDP_prec
  )

  if (plot_graph) {
    plot <- ggplot(df, aes(x = time_f, y = MDP_mean)) +
      geom_errorbar(aes(ymin = MDP_low, ymax = MDP_high), alpha = .1, color = "red") +
      geom_point(size = .1)+
      #scale_x_continuous(limits = c(-9000, 2000), breaks = seq(-9000, 2000, by = 2000), minor_breaks = seq(-4000, -3000, by = 10)) +
      scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,1), minor_breaks = seq(-1,11, by = 0.1))+
      theme_pubr() +
      labs(
        title = "Marginal Data Posterior Output",
        x = "Time AD",
        y = "Mean MDPs"
      )
    print(plot)
  }

  if (print) {
    print(df)
  }
  class(MDPs) <- ("MDPs")
  return(MDPs)
}
