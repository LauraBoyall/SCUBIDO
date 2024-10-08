SCUBIDO_cal.sorted <- function(sorted, plot = TRUE, summary = TRUE) {
  # Function body remains unchanged here
  cal_data <- list(
    N.rows_m = nrow(sorted$xrf_m),
    N.cols_m = ncol(sorted$xrf_m),
    temp_m = sorted$temp_m,
    xrf_m = sorted$xrf_m,
    time_m = sorted$time_m,
    k = ncol(sorted$xrf_m),
    R = diag(ncol(sorted$xrf_m))
  )

  cal_pars_to_save <- c(
    "beta0", "beta1", "beta2",
    "Sigma_inv",
    "xrf_m_pred",
    "sd_rw",
    "mean_m"
  )

  jags_file <- system.file("jags_models", "modern_MVN_polynomial_model.jags", package = "SCUBIDO")

  cal_file <- paste0("cal_model_", Sys.Date(), ".rds")
  if (!file.exists(cal_file)) {
    cal_run <- jags(
      data = cal_data,
      parameters.to.save = cal_pars_to_save,
      model.file = jags_file
    )
    saveRDS(cal_run, file = cal_file)
  } else {
    cal_run <- readRDS(file = cal_file)
  }

  calibration <- list (
    sims.list = cal_run$BUGSoutput$sims.list,
    sorted = sorted
  )

  if (plot) {
    par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
    par(mfrow = c(7,2))
    par_means = cal_run$BUGSoutput$mean
    temp_grid = seq(-3,3,length=50)
    for(i in 1:11) {
      plot(cal_data$temp_m, cal_data$xrf_m[,i],
           pch = 19,
           xlab = 'Temperature',
           ylab = 'XRF',
           main = colnames(cal_data$xrf_m)[i],
           ylim = range(c(cal_data$xrf_m[,i], cal_data$xrf_m[,i])))

      mean_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2))
      lines(temp_grid, mean_pred)

      low_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2) - 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
      lines(temp_grid, low_pred, lty = 2)

      high_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2) + 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
      lines(temp_grid, high_pred, lty = 2)
    }
    par(mfrow = c(1,1))
  }

  if (summary) {
    print(cal_run)
  }

  class(calibration) <- ("calibration")
  return(calibration)
}
