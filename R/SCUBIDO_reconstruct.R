#' Reconstruction of climate
#'
#' SCUBIDO_reconstruct is a function which sorts the previous analysis from the
#'     calibration period and then applies what is learnt on the fossil_data.
#'     This first creates marginal data posteriors to transform the multivariate
#'     proxy data into a single estimation of climate based on each time point
#'     of the XRF data. A random walk is then applied to join each of these layers
#'     and reconstruct climate through time and quantifies 95% and 50% confidence
#'     intervals.
#'     Note that this function can take many hours to run depending on how many data
#'     points in the XRF data, a dataset containing approximately 60,000 time points
#'     takes around 12 hours on a standard computer.
#'
#' @import ggplot2
#' @import zoo
#' @import R2jags
#' @import mvtnorm
#' @import fitdistrplus
#' @import ggpubr
#'
#' @param calibration_data is the output list from the SCUBIDO_cal function.
#' @param plot_graph A flag to plot the output of the reconstruction.
#'
#' @return a full reconstruction of climate given the proxy data.
#'
#' @examples
#' \dontrun{
#' SCUBIDO_reconstruct(calibration_data, plot_graph = TRUE)
#' }
#'

SCUBIDO_reconstruct <- function(calibration_data,  plot_graph = FALSE) {

  # Display a message warning the user about the potential long runtime
  message("Warning: This operation may take a significant amount of time to complete.")

  # Ask the user for confirmation to proceed
  user_input <- readline(prompt = "Press [Enter] to continue or [Esc] to cancel the operation: ")

  temp_grid <- seq(-3,3, length = 50)

  # Part 1: SCUBIDO_apply functionality
  MDP_file <- paste0("MDPs_mising", Sys.Date(), ".rds")

  cal_run_results <- calibration_data$sims.list
  N_rows_f <- nrow(calibration_data$sorted$xrf_f_resc)

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
            dmvnorm(calibration_data$sorted$xrf_f_resc[i, ],
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
    time_f = calibration_data$sorted$time_f,
    MDP_mean = MDP_mean,
    MDP_low = MDP_low,
    MDP_high = MDP_high,
    MDP_prec = MDP_prec,
    sd_rw = cal_run_results,
    sorted = calibration_data$sorted,
    calibration = calibration_data$sims.list
  )

  # Part 2: SCUBIDO_reconstruct functionality
  time_fm <- MDPs$sorted$time_fm
  time_all <- MDPs$sorted$time_all
  time_grid <- MDPs$time_grid

  pick <- which(time_all %in% time_fm)

  MDP_mean_all <- c(MDPs$MDP_mean[,1], MDPs$sorted$temp_m)
  MDP_prec_all <- c(MDPs$MDP_prec[, 1], rep(1e5, length(MDPs$sorted$temp_m)))

  lnorm_pars <- fitdist(as.vector(MDPs$calibration$sd_rw), "lnorm")

  fossil_model_data <- list(
    N.rows_fm = length(MDP_mean_all),
    N.rows_all = length(time_all),
    pick = pick,
    MDP_mean = MDP_mean_all,
    MDP_prec = MDP_prec_all,
    time_all = time_all,
    a_rw = lnorm_pars$estimate[1],
    b_rw = 1 / (lnorm_pars$estimate[2]^2)
  )

  fossil_pars_to_watch <- c("temp_all", "sd_rw")
  jags_file <- system.file("jags_models", "fossil_MDP_model.jags", package = "SCUBIDO")

  fossil_file <- paste0("fossil_model_missing", Sys.Date(), ".rds")
  if (!file.exists(fossil_file)) {
    fossil_run <- jags(
      data = fossil_model_data,
      parameters.to.save = fossil_pars_to_watch,
      model.file = jags_file,
      n.chains = 4,
      n.iter = 6000,
      n.burnin = 4000,
      n.thin = 4
    )
    saveRDS(fossil_run, fossil_file)
  } else {
    fossil_run <- readRDS(file = fossil_file)
  }

  temp_all <- fossil_run$BUGSoutput$sims.list$temp_all
  pick_grid <- which(time_all %in% MDPs$sorted$time_grid)
  temp_grid <- temp_all[, pick_grid]

  df_final <- data.frame(
    Age_BP = MDPs$sorted$time_grid,
    Age_AD =  1950 - MDPs$sorted$time_grid,
    temp_med = apply(temp_grid, 2, median),
    temp_low_95 = apply(temp_grid, 2, quantile, 0.025),
    temp_high_95 = apply(temp_grid, 2, quantile, 0.975),
    temp_low_50 = apply(temp_grid, 2, quantile, 0.25),
    temp_high_50 = apply(temp_grid, 2, quantile, 0.75),
    temp_low_75 = apply(temp_grid, 2, quantile, 0.125),
    temp_high_75 = apply(temp_grid, 2, quantile, 0.875)
  )

  if (plot_graph) {
    smot <- 1
    plot <- ggplot(df_final, aes(x = yBP / 1000)) +
      geom_ribbon(aes(ymin = rollmean(temp_low_95, smot, na.pad = TRUE), ymax = rollmean(temp_high_95, smot, na.pad = TRUE)), alpha = 0.6, fill = 'lightblue') +
      geom_ribbon(aes(ymin = rollmean(temp_low_50, smot, na.pad = TRUE), ymax = rollmean(temp_high_50, smot, na.pad = TRUE)), alpha = 0.3, fill = 'darkblue') +
      geom_line(aes(y = rollmean(temp_med, smot, na.pad = TRUE)), linetype = 'solid', size = 0.2) +
      geom_line(aes(y = rollmean(temp_high_95, smot, na.pad = TRUE)), linetype = 'dashed', size = 0.05) +
      geom_line(aes(y = rollmean(temp_low_95, smot, na.pad = TRUE)), linetype = 'dashed', size = 0.05) +
      geom_hline(yintercept = 0) + ggpubr::theme_pubr() +
      labs(
        title = "Palaeoclimate Reconstruction",
        x = "Age (years BP)",
        y = expression("Temperature "(degree * C))
      )
    print(plot)
  }

  climate_reconstruction <- list(
    reconstruction = df_final,
    model_output = fossil_run
  )

  return(climate_reconstruction)
}

