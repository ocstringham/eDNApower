#' Compute Frequentist Power for eDNA Samples
#'
#' Computes the probability of observing at least `s_min` positive samples
#' given the number of samples, number of technical replicates, and detection probabilities.
#'
#' @param n Integer. Number of samples collected.
#' @param k Integer. Number of technical replicates per sample.
#' @param theta Numeric between 0 and 1. Probability that the species is present in a sample.
#' @param p Numeric between 0 and 1. Probability that a technical replicate detects the species if it is present.
#' @param s_min Integer. Minimum number of positive samples to be considered successful (default 5).
#'
#' @return Numeric. The frequentist power (probability of ≥ s_min positive samples).
#' @examples
#' freq_power(n = 10, k = 3, theta = 0.6, p = 0.8, s_min = 5)
#' freq_power(n = 15, k = 2, theta = 0.7, p = 0.9)
#' @export
freq_power <- function(n, k, theta, p, s_min = 5) {
  pi <- theta * (1 - (1 - p)^k)
  1 - pbinom(s_min - 1, size = n, prob = pi)
}

#' Compute Frequentist Power Grid
#'
#' Computes frequentist power across a grid of sample sizes (`n_vals`) and technical replicates (`k_vals`).
#'
#' @param n_vals Integer vector. Possible number of samples.
#' @param k_vals Integer vector. Possible numbers of technical replicates per sample.
#' @param theta Numeric between 0 and 1. Probability species is present in a sample.
#' @param p Numeric between 0 and 1. Probability a technical replicate detects the species if present.
#' @param s_min Integer. Minimum number of positive samples to be considered successful (default 5).
#'
#' @return Data frame with columns `n`, `k`, and `power`.
#' @examples
#' freq_power_grid(5:10, 1:3, theta = 0.6, p = 0.8)
#' @export
freq_power_grid <- function(n_vals, k_vals, theta, p, s_min = 5) {
  grid <- expand.grid(n = n_vals, k = k_vals)
  grid$power <- mapply(freq_power, n = grid$n, k = grid$k,
                       MoreArgs = list(theta = theta, p = p, s_min = s_min))
  return(grid)
}

#' Plot Frequentist Power Heatmap
#'
#' Generates a heatmap of frequentist power across combinations of `n` and `k`.
#'
#' @param grid Data frame returned by `freq_power_grid`.
#'
#' @return ggplot2 object displaying a heatmap.
#' @examples
#' grid <- freq_power_grid(5:10, 1:3, theta = 0.6, p = 0.8)
#' plot_freq_heatmap(grid)
#' @export
plot_freq_heatmap <- function(grid) {
  ggplot2::ggplot(grid, ggplot2::aes(x = n, y = k, fill = power)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_fermenter(limits = c(0, 1), breaks = c(0.8,0.9,0.95, 0.99)) +
    ggplot2::labs(x = "Number of samples (n)",
         y = "Technical replicates per sample (k)",
         fill = "Frequentist Power",
         title = "Frequentist Power Heatmap") +
    ggplot2::theme_minimal()
}
