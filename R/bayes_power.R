#' Compute Bayesian Prior Predictive Power for eDNA Samples
#'
#' Simulates the expected probability of observing at least `s_min` positive samples
#' accounting for uncertainty in sample availability (`theta`) and replicate detection (`p`).
#'
#' @param n Integer. Number of samples collected.
#' @param k Integer. Number of technical replicates per sample.
#' @param s_min Integer. Minimum number of positive samples to be considered successful (default 5).
#' @param n_sim Integer. Number of Monte Carlo simulations to run (default 10000).
#' @param theta_prior Numeric vector. Optional posterior or prior draws for theta. Either this or a_theta/b_theta must be provided.
#' @param p_prior Numeric vector. Optional posterior or prior draws for p. Either this or a_p/b_p must be provided.
#' @param a_theta Numeric. Alpha parameter for Beta prior of theta (optional).
#' @param b_theta Numeric. Beta parameter for Beta prior of theta (optional).
#' @param a_p Numeric. Alpha parameter for Beta prior of p (optional).
#' @param b_p Numeric. Beta parameter for Beta prior of p (optional).
#'
#' @return Numeric. Estimated Bayesian prior predictive power (probability of ≥ s_min positive samples).
#' @examples
#' bayes_power(n = 10, k = 3, s_min = 5, a_theta = 6, b_theta = 4, a_p = 8, b_p = 2)
#' theta_post <- rbeta(10000, 6, 4)
#' p_post <- rbeta(10000, 8, 2)
#' bayes_power(n = 10, k = 3, s_min = 5, theta_prior = theta_post, p_prior = p_post)
#' @export
bayes_power <- function(n, k, s_min = 5, n_sim = 10000,
                        theta_prior = NULL, p_prior = NULL,
                        a_theta = NULL, b_theta = NULL,
                        a_p = NULL, b_p = NULL) {

  # theta draws
  if (!is.null(theta_prior)) {
    theta_sim <- theta_prior
    n_sim <- length(theta_sim)
  } else if (!is.null(a_theta) & !is.null(b_theta)) {
    theta_sim <- rbeta(n_sim, a_theta, b_theta)
  } else stop("Provide theta_prior vector or Beta params a_theta/b_theta")

  # p draws
  if (!is.null(p_prior)) {
    p_sim <- p_prior
    n_sim <- length(p_sim)
  } else if (!is.null(a_p) & !is.null(b_p)) {
    p_sim <- rbeta(n_sim, a_p, b_p)
  } else stop("Provide p_prior vector or Beta params a_p/b_p")

  pi_sim <- theta_sim * (1 - (1 - p_sim)^k)
  S_sim <- rbinom(n_sim, size = n, prob = pi_sim)
  mean(S_sim >= s_min)
}

#' Compute Bayesian Power Grid
#'
#' Computes Bayesian prior predictive power across a grid of sample sizes (`n_vals`) and technical replicates (`k_vals`).
#'
#' @param n_vals Integer vector. Possible number of samples.
#' @param k_vals Integer vector. Possible numbers of technical replicates per sample.
#' @param s_min Integer. Minimum number of positive samples to be considered successful (default 5).
#' @param n_sim Integer. Number of Monte Carlo simulations to run (default 10000).
#' @param theta_prior Numeric vector. Optional posterior or prior draws for theta.
#' @param p_prior Numeric vector. Optional posterior or prior draws for p.
#' @param a_theta Numeric. Alpha parameter for Beta prior of theta (optional).
#' @param b_theta Numeric. Beta parameter for Beta prior of theta (optional).
#' @param a_p Numeric. Alpha parameter for Beta prior of p (optional).
#' @param b_p Numeric. Beta parameter for Beta prior of p (optional).
#'
#' @return Data frame with columns `n`, `k`, and `power`.
#' @examples
#' bayes_power_grid(5:10, 1:3, s_min = 5, a_theta = 6, b_theta = 4, a_p = 8, b_p = 2)
#' @export
bayes_power_grid <- function(n_vals, k_vals, s_min = 5, n_sim = 10000,
                             theta_prior = NULL, p_prior = NULL,
                             a_theta = NULL, b_theta = NULL,
                             a_p = NULL, b_p = NULL) {
  grid <- expand.grid(n = n_vals, k = k_vals)
  grid$power <- mapply(bayes_power, n = grid$n, k = grid$k,
                       MoreArgs = list(s_min = s_min, n_sim = n_sim,
                                       theta_prior = theta_prior, p_prior = p_prior,
                                       a_theta = a_theta, b_theta = b_theta,
                                       a_p = a_p, b_p = b_p))
  return(grid)
}

#' Plot Bayesian Power Heatmap
#'
#' Generates a heatmap of Bayesian prior predictive power across combinations of `n` and `k`.
#'
#' @param grid Data frame returned by `bayes_power_grid`.
#'
#' @return ggplot2 object displaying a heatmap.
#' @examples
#' grid <- bayes_power_grid(5:10, 1:3, s_min = 5, a_theta = 6, b_theta = 4, a_p = 8, b_p = 2)
#' plot_bayes_heatmap(grid)
#' @export
plot_bayes_heatmap <- function(grid) {
  ggplot2::ggplot(grid, ggplot2::aes(x = n, y = k, fill = power)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_fermenter(limits = c(0, 1), breaks = c(0.8,0.9,0.95, 0.99)) +
    ggplot2::labs(x = "Number of samples (n)",
         y = "Technical replicates per sample (k)",
         fill = "Bayesian Power",
         title = "Bayesian Prior Predictive Power Heatmap") +
    ggplot2::theme_minimal()
}
