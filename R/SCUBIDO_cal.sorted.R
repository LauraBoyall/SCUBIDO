#' Calculate the modern relationship between the proxy and climate
#'
#' The `SCUBIDO_cal` function uses the modern data from the \code{\link{SCUBIDO_input}}
#'     function and computes the relationship between the modern XRF data and an
#'     instrumental climate time series. The function uses the JAGS package to fit
#'     a multivariate polynomial regression model to identify the relationship
#'     between the different XRF elements and the climate data provided. Specific
#'     model details can be found in Boyall et al (in prep). The relationship between
#'     the elements and climate is then expressed through a likelihood function.
#'     If check_convergence = TRUE a plot showing the Rhat values is produced. This
#'     checks whether the Marcov Chain Monte Carlo (MCMC) algorithm has fitted.
#'     If a point is <1.05 then is is assumed that the model has converged well.
#'     The saved results from this function will be used to form the final reconstruction
#'     in the \code{\link{SCUBIDO_reconstruct}} function.
#'
#'
#' @param sorted the modern dataset saved after using the \code{\link{SCUBIDO_input}} function
#' @param plot returns a plot of the relationship between the modern XRF elements and climate
#' @param summary returns a printed summary of the output of the calibration model
#' @export
#'
#' @import R2jags
#'
#' @examples
#' \dontrun{
#' SCUBIDO_cal(x,  plot = TRUE, summary = TRUE)
#' }
#'

SCUBIDO_cal<- function(sorted,  plot = TRUE, summary = TRUE){
  UseMethod("SCUBIDO_cal")
}
#' @export
SCUBIDO_cal.sorted <- function(sorted,  plot = TRUE, summary = TRUE) {
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
           ylim = range(c(cal_data$xrf_m[,i], cal_data$xrf_m[,i]))) #from the for loop to the end of this, it is just plotting the NAO against the XRF, nothing to do with the model

      mean_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2))
      lines(temp_grid, mean_pred) #Here we are putting the parameter means on each of the graphs.

      low_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2) - 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
      lines(temp_grid, low_pred, lty = 2)

      high_pred = with(par_means, beta0[i] + beta1[i] * temp_grid + beta2[i] * (temp_grid^2) + 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
      lines(temp_grid, high_pred, lty = 2)
      #abline(h = xrf_f_resc[510,i], col = 'red')
    }
    par(mfrow = c(1,1))
  }


  if (summary) {
    print(cal_run)
  }

  class(calibration) <- ("calibration")
  return(calibration)
}

