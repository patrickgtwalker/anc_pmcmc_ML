plot_eir <- function(history, true_history, times, obs_end = NULL) {
  if (is.null(obs_end)) {
    obs_end <- max(times)
  }
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  matplot(times, t(history['EIR', , -1]), type = "l",
          xlab = "Time", ylab = "EIR",
          col = "#A6CEE3", lty = 1, ylim = range(history['EIR', , -1]))
  matpoints(times, true_history$EIR_true, pch = 19,
            col = "#1F78B4")
}
