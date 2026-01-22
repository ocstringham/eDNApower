#' Compute Power Analysis for eDNA Hierarchical Detection Sampling Design
#'
#' Computes the probability of observing at least `s_min` positive samples
#' given the number of samples, number of technical replicates, and values or distributions for detection probabilities.
#' This function considers at least one technical replicate per sample detecting DNA as a positive detection.
#'
#' @param n Integer. Number of samples collected.
#' @param k Integer. Number of technical replicates per sample.
#' @param theta Numeric or Numeric Vector. Value(s) between 0 and 1. Probability that DNA is present in a sample.
#' @param p Numeric or Numeric Vector. Value(s) between 0 and 1. Probability that a technical replicate detects the DNA if it is present.
#' @param s_min Integer. Minimum number of positive samples to be considered successful.
#' @param n_sim Integer. Number of Monte Carlo simulations to run (default 10,000). Only used if a vector is provided for `theta` or `p`.
#'
#' @return Numeric. The probability of ≥ s_min positive samples, considering at least one technical replicate per sample detecting DNA as a positive detection.
#' @examples
#' eDNA_power(n = 10, k = 3, theta = 0.6, p = 0.8, s_min = 5)
#' eDNA_power(n = 15, k = 2, theta = runif(1000, 0.4, 0.6), p = rbeta(1000, 1.8, 4.2), s_min = 5)
#' @export
eDNA_power <- function(n, k, theta, p, s_min, n_sim = 10000) {

  # if theta a p are single values
  if (length(theta) == 1 && length(p) == 1) {
    pi <- theta * (1 - (1 - p)^k)
    1 - pbinom(s_min - 1, size = n, prob = pi)

  # if theta or p are vectors, run Monte Carlo simulations
  } else{
    pi_sim <- theta * (1 - (1 - p)^k)
    S_sim <- rbinom(n_sim, size = n, prob = pi_sim)
    mean(S_sim >= s_min)
  }
}



#' Compute Power Grid
#'
#' Computes prior predictive power across a grid of sample sizes (`n_vals`) and technical replicates (`k_vals`).
#'
#' @param n_vals Integer vector. Possible number of samples.
#' @param k_vals Integer vector. Possible numbers of technical replicates per sample.
#' @param theta Numeric or Numeric Vector. Value(s) between 0 and 1. Probability that DNA is present in a sample.
#' @param p Numeric or  Numeric Vector. Value(s) between 0 and 1. Probability that a technical replicate detects the DNA if it is present.
#' @param s_min Integer. Minimum number of positive samples to be considered successful.
#' @param n_sim Integer. Number of Monte Carlo simulations to run (default 10,000). Only used if a vector is provided for `theta` or `p`.
#'
#' @return Data frame with columns `n`, `k`, and `power`.
#' @examples
#' eDNA_power_grid(n_vals = 5:10, k_vals = 1:3, theta = 0.6, p = 0.8, s_min = 5)
#' eDNA_power_grid(n_vals = 5:10, k_vals = 1:3, theta = runif(1000, 0.4, 0.6), p = rbeta(1000, 1.8, 4.2), s_min = 5)
#' @export
eDNA_power_grid <- function(n_vals, k_vals, s_min, n_sim = 10000,
                             theta = NULL, p = NULL) {
  grid <- expand.grid(n = n_vals, k = k_vals)
  grid$power <- mapply(eDNA_power, n = grid$n, k = grid$k,
                       MoreArgs = list(s_min = s_min, n_sim = n_sim,
                                       theta = theta, p = p))
  return(grid)
}



#' Extract n sample and k replicates that acheive target power
#'
#' Extracts the minimum number of samples (n) and technical replicates (k) needed to achieve a target power level.
#'
#' @param grid Data frame returned by `eDNA_power_grid`.
#' @param target_power Numeric. Target power level between 0 and 1.
#' @return Data frame with minimum `n` and `k` to achieve target power.
#' @examples
#' grid <- eDNA_power_grid(n_vals = 5:20, k_vals = 1:5, theta = 0.6, p = 0.8, s_min = 5)
#' eDNA_power_extract(grid)
#'
#'
#' @export
eDNA_power_extract <- function(grid, target_power){


}


#' Plot Power Heatmap
#'
#' Generates a heatmap of prior predictive power across combinations of `n` and `k`.
#'
#' @param grid Data frame returned by `eDNA_power_grid`.
#' @param breaks Numeric vector. Break points for the color (probability) scale (default c(0.5, 0.75, 0.95, 0.99)).
#' @return ggplot2 object displaying a heatmap.
#' @examples
#' grid <- eDNA_power_grid(n_vals = 5:10, k_vals = 1:3, theta = 0.6, p = 0.8, s_min = 5)
#' eDNA_power_heatmap(grid)
#' @export
eDNA_power_heatmap <- function(grid, breaks = c(0.5,0.75,0.95, 0.99)) {
  ggplot2::ggplot(grid, ggplot2::aes(x = n, y = k, fill = power)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_fermenter(limits = c(0, 1), breaks = breaks, palette = "YlGn", direction = 1) +
    ggplot2::labs(x = "Number of samples (n)",
                  y = "Technical replicates per sample (k)",
                  fill = "Probability of ≥ s_min\npositive samples") +
    ggplot2::theme_minimal()
}
