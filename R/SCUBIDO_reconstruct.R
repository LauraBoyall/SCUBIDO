#' Reconstruction of climate
#'
#' SCUBIDO_reconstruct is a function which sorts the previous analysis from the
#'     calibration period and the marginal data posteriors and applies a random
#'     walk
#'
#' @import ggplot2
#' @import zoo
#' @import R2jags
#'
#' @param MDPs the MDP list created in the apply function
#' @param time_grid The maximum range of age in BCE/CE from your modern data set to your fossil data set e.g. range(fossil_data$age, modern_data$age) and then sequence by wanted resolution e.g. 10 years. see example.
#' @param plot A plot of the final reconstruction
#'
#' @return a full reconstruction of climate given the proxt data
#'
#'
#' @examples
#' \dontrun{
#' time_grid <- seq(-8331.625, 1932.000, by = 10)
#' SCUBIDO_reconstruct(time_grid, plot = TRUE)
#' }
#'
SCUBIDO_reconstruct.MDPs <- function(MDPs, time_grid, plot_graph = FALSE){
UseMethod("SCUBIDO_reconstruct")
}
#'
#'
#'
#'
#' @export
#'
SCUBIDO_reconstruct <- function(MDPs,time_grid, plot = FALSE) {
  # Put all the times together
  time_fm <- sort(c(MDPs$sorted$time_f, MDPs$sorted$time_m))
  time_all <- sort(c(time_fm, time_grid))

  # Avoid having any repeating times
  which(diff(time_all) == 0)  # Should be integer(0)

  # Find the times in the grid which have data
  pick <- which(time_all %in% time_fm)

  # Set up matrix
  MDP_mean_all <- c(MDPs$MDP_mean[,1], MDPs$sorted$temp_m)

  MDP_prec_all <- c(MDPs$MDP_prec[, 1], rep(1e5, length(MDPs$sorted$temp_m)))

  # Get prior estimates for the random walk varianc
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
  # Check the fit
  # plot(fossil_run)

  # Now plot the overall run
  temp_all <- fossil_run$BUGSoutput$sims.list$temp_all

  pick_grid <- which(time_all %in% time_grid)
  temp_grid <- temp_all[, pick_grid]


  df_final <- data.frame(
    time_grid = time_grid,
    temp_med = apply(temp_grid, 2, median),
    temp_low_95 = apply(temp_grid, 2, quantile, 0.025),
    temp_high_95 = apply(temp_grid, 2, quantile, 0.975),
    temp_low_50 = apply(temp_grid, 2, quantile, 0.25),
    temp_high_50 = apply(temp_grid, 2, quantile, 0.75)
  )

  df_final$yBP <- 1950 - df_final$time_grid
  df_final$yAD <- df_final$time_grid

  if(plot){
   # library(zoo)
    #library(ggplot2)

    smot <- 1 #rolling mean smoothing level

    plot <- ggplot(df_final, aes(x =yBP/1000))+
      geom_ribbon(aes(ymin = rollmean(`temp_low_95`, smot, na.pad = TRUE), ymax = rollmean(`temp_high_95`, smot, na.pad = TRUE)), alpha = 0.6, fill = 'lightblue') +
      geom_ribbon(aes(ymin = rollmean(`temp_low_50`, smot, na.pad = TRUE), ymax = rollmean(`temp_high_50`, smot, na.pad = TRUE)), alpha = 0.3, fill = 'darkblue') +
      geom_line(aes(y = rollmean(temp_med, smot, na.pad = TRUE)), linetype = 'solid', size = 0.2) +
      geom_line(aes(y = rollmean(`temp_high_95`, smot, na.pad = TRUE)),  linetype = 'dashed',size = 0.05) +
      geom_line(aes(y = rollmean(`temp_low_95`, smot, na.pad = TRUE)), linetype = 'dashed',  size = 0.05) +
      geom_hline(yintercept = 0) + ggpubr::theme_pubr() +
      #scale_x_reverse(limits = c(10.5,0), breaks = seq(10.5,0, -2.5), minor_breaks = seq(4,5, by = 0.05)) +
      scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,.5), minor_breaks = seq(-1,11, by = 0.1))+
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

  # Assign the result_list to a variable in the global environment
  #assign("climate_reconstruction", climate_reconstruction, envir = .GlobalEnv)
  return(climate_reconstruction)


}


#final <- reconstruction(time_grid, plot = T)


save_all<- function(modern_data, fossil_data){

  save(modern_data, fossil_data, result_list, MDPs, calibration_list, climate_reconstruction,
       file = "Climate reconstruction.RData")

}
