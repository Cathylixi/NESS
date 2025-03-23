#' Plot Stability or Rareness Scores (Improved)
#'
#' Generates a log-scaled line plot with improved aesthetics for visualizing metrics like stability or rareness score over GCP.
#'
#' @param x Numeric vector of GCP values (x-axis).
#' @param y Numeric vector of corresponding metric values (y-axis).
#' @param ylab Character label for the y-axis.
#' @param main_title Character string for the plot title.
#'
#' @return A base R plot.
#' @keywords internal
plot_improved <- function(x, y, ylab, main_title) {
  plot(x, y,
       log = "x",
       main = main_title,
       xlab = "GCP",
       ylab = ylab,
       type = "l",
       col = "black",
       lwd = 2,
       cex.lab = 1.4,
       cex.main = 1.6,
       cex.axis = 1.2
  )
  grid(col = "gray", lty = "dotted")
  box(lwd = 2)
}
