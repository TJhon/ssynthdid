# estimation.R
# Main estimation functions for staggered synthetic control

#' Compute treatment effect estimate
#'
#' @param Y_adjusted Outcome matrix adjusted for covariates
#' @param omega Unit weights
#' @param lambda Time weights
#' @param N1 Number of treated units
#' @param T1 Number of post-treatment periods
#' @return Scalar treatment effect estimate
compute_treatment_effect <- function(Y_adjusted, omega, lambda, N1, T1) {
  omega_full <- c(-omega, rep(1 / N1, N1))
  lambda_full <- c(-lambda, rep(1 / T1, T1))

  as.numeric(crossprod(omega_full, Y_adjusted %*% lambda_full))
}
