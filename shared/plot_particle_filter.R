plot_particle_filter <- function(history, true_history, times) {
  cis <- addCIs(true_history,true_history$positive,true_history$tested)
  cols <- ncol(cis)
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  matplot(times, t(history[1, , -1]), type = "l",
          xlab = "Time", ylab = "Prevalence",
          col = "#A6CEE3", lty = 1, ylim = range(c(cis[,(cols-2):cols],history[1, , -1]),na.rm=TRUE))
  matpoints(times, cis$mean, pch = 19,
            col = "#1F78B4")
  arrows(cis$t, cis$lower, cis$t, cis$upper, length=0.05, angle=90, code=3, col = "#1F78B4")
}
