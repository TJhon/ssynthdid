#' Contract list of matrices with coefficient vector
#'
#' @param X List of matrices with same dimensions
#' @param v Coefficient vector
#' @return Weighted sum of matrices

contract_matrices <- function(X, v) {
# contract3 <- function(X, v) {
  # same size
  X1 <- X[[1]]
  n <- nrow(X1)
  t <- ncol(X1)
  out <- matrix(0, nrow = n, ncol = t)
  if (length(v) == 0) { return(out) }

  for (ii in 1:length(v)) {
    out = out + v[ii] * X[[ii]]
  }
  return(out)
}


#' Frank-Wolfe step for regularized least squares on simplex
#'
#' @param A Design matrix
#' @param x Current weight vector (on simplex)
#' @param b Target vector
#' @param eta Ridge penalty parameter
#' @param alpha Fixed step size (if NULL, uses line search)
#' @return Updated weight vector
fw_step <- function(A, x, b, eta, alpha = NULL) {

  # Compute gradient at current point
  Ax <- A %*% x
  half_grad <- crossprod(Ax - b, A) + eta * x

  # Find extreme point of simplex in descent direction
  i <- which.min(half_grad)

  # Fixed step size update
  if (!is.null(alpha)) {
    x <- x * (1 - alpha)
    x[i] <- x[i] + alpha
    return(x)
  }

  # Line search update
  direction <- -x
  direction[i] <- 1 - x[i]

  # Check if already at optimum
  if (all(direction == 0)) {
    return(x)
  }

  # Compute optimal step size
  residual_direction <- A[, i, drop = FALSE] - Ax
  numerator <- -as.numeric(crossprod(half_grad, direction))
  denominator <- sum(residual_direction^2) + eta * sum(direction^2)

  step_size <- numerator / denominator
  constrained_step <- min(1, max(0, step_size))

  x + constrained_step * direction
}


#' Solve for synthetic control weights using Frank-Wolfe
#'
#' @param Y Outcome matrix (N0 x T) where last column is target
#' @param zeta Regularization parameter (complex for consistency)
#' @param intercept Whether to demean outcomes
#' @param lambda Initial weights (if NULL, uniform)
#' @param min_decrease Convergence threshold
#' @param max_iter Maximum iterations
#' @return List with optimized weights and objective values
sc_weight_fw <- function(Y, zeta, intercept = TRUE, lambda = NULL,
                         min_decrease = 1e-3, max_iter = 1000) {

  T0 <- ncol(Y) - 1
  N0 <- nrow(Y)

  # Validate inputs
  if (T0 <= 0) {
    stop(sprintf("Invalid Y matrix: need at least 2 columns, got %d", ncol(Y)))
  }
  if (N0 <= 0) {
    stop(sprintf("Invalid Y matrix: need at least 1 row, got %d", N0))
  }

  # Initialize weights
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }

  # Demean if requested
  if (intercept) {
    Y <- apply(Y, 2, function(col) col - mean(col))
  }

  # Setup optimization
  A <- Y[, seq_len(T0), drop = FALSE]
  b <- Y[, T0 + 1]
  eta <- N0 * Re(zeta^2)

  # Optimization loop
  iter <- 0
  objective_vals <- numeric(0)

  while (iter < max_iter && (iter < 2 || objective_vals[iter - 1] - objective_vals) ) {
    iter <- iter + 1

    # Frank-Wolfe step
    lambda <- fw_step(A, lambda, b, eta)

    # Compute objective
    residual <- Y[seq_len(N0), ] %*% c(lambda, -1)
    obj_val <- Re(zeta^2) * sum(lambda^2) + sum(residual^2) / N0
    objective_vals <- c(objective_vals, obj_val)

    # Check convergence
    if (iter >= 2) {
      improvement <- objective_vals[iter - 1] - objective_vals[iter]
      if (improvement < min_decrease^2) {
        break
      }
    }
  }

  list(lambda = lambda, vals = objective_vals)
}


#' Solve for synthetic control with covariates using Frank-Wolfe
#'
#' @param Y Outcome matrix (collapsed form)
#' @param X List of covariate matrices (collapsed form)
#' @param zeta_lambda Regularization for time weights
#' @param zeta_omega Regularization for unit weights
#' @param lambda_intercept Whether to demean for lambda optimization
#' @param omega_intercept Whether to demean for omega optimization
#' @param min_decrease Convergence threshold
#' @param max_iter Maximum iterations
#' @param lambda Initial time weights
#' @param omega Initial unit weights
#' @param beta Initial covariate coefficients
#' @param update_lambda Whether to update lambda
#' @param update_omega Whether to update omega
#' @return List with optimized weights and coefficients
sc_weight_fw_covariates <- function(Y, X = list(),
                                    zeta_lambda = 0, zeta_omega = 0,
                                    lambda_intercept = TRUE,
                                    omega_intercept = TRUE,
                                    min_decrease = 1e-3,
                                    max_iter = 1000,
                                    lambda = NULL,
                                    omega = NULL,
                                    beta = NULL,
                                    update_lambda = TRUE,
                                    update_omega = TRUE) {

  T0 <- ncol(Y) - 1
  N0 <- nrow(Y) - 1

  # Initialize parameters
  if (is.null(lambda)) lambda <- rep(1 / T0, T0)
  if (is.null(omega)) omega <- rep(1 / N0, N0)
  if (is.null(beta)) beta <- rep(0, length(X))

  # Function to update weights given adjusted outcomes
  update_weights <- function(Y_adjusted, lambda_current, omega_current) {

    # Validate dimensions
    if (N0 <= 0 || T0 <= 0) {
      stop("Invalid dimensions in update_weights")
    }

    # Update lambda (time weights)
    Y_lambda_full <- Y_adjusted[seq_len(N0), , drop = FALSE]
    Y_lambda <- if (lambda_intercept) {
      apply(Y_lambda_full, 2, function(col) col - mean(col))
    } else {
      Y_lambda_full
    }

    if (update_lambda) {
      lambda_current <- fw_step(
        Y_lambda[, seq_len(T0), drop = FALSE],
        lambda_current,
        Y_lambda[, T0 + 1],
        N0 * Re(zeta_lambda^2)
      )
    }

    residual_lambda <- Y_lambda %*% c(lambda_current, -1)

    # Update omega (unit weights)
    Y_omega_full <- t(Y_adjusted[, seq_len(T0), drop = FALSE])
    Y_omega <- if (omega_intercept) {
      apply(Y_omega_full, 2, function(col) col - mean(col))
    } else {
      Y_omega_full
    }

    if (update_omega) {
      omega_current <- fw_step(
        Y_omega[, seq_len(N0), drop = FALSE],
        omega_current,
        Y_omega[, N0 + 1],
        T0 * Re(zeta_omega^2)
      )
    }

    residual_omega <- Y_omega %*% c(omega_current, -1)

    # Compute objective value
    obj_val <- Re(zeta_omega^2) * sum(omega_current^2) +
      Re(zeta_lambda^2) * sum(lambda_current^2) +
      sum(residual_omega^2) / T0 +
      sum(residual_lambda^2) / N0

    list(
      val = obj_val,
      lambda = lambda_current,
      omega = omega_current,
      residual_lambda = residual_lambda,
      residual_omega = residual_omega
    )
  }

  # Gradient function for beta (using environment for efficiency)
  compute_beta_gradient <- function(Xi, weights) {
    # Safe indexing for both dimensions
    Xi_control <- Xi[seq_len(N0), , drop = FALSE]
    Xi_pre <- Xi[, seq_len(T0), drop = FALSE]

    as.numeric(
      crossprod(weights$residual_lambda, Xi_control %*% c(weights$lambda, -1)) / N0 +
        crossprod(weights$residual_omega, t(Xi_pre) %*% c(weights$omega, -1)) / T0
    )
  }

  # Main optimization loop
  iter <- 0
  objective_vals <- numeric(0)

  Y_adjusted <- Y - contract_matrices(X, beta)
  weights <- update_weights(Y_adjusted, lambda, omega)

  while (iter < max_iter) {
    iter <- iter + 1

    # Update beta using gradient descent with diminishing step size
    if (length(X) > 0) {
      grad_beta <- vapply(X, compute_beta_gradient, numeric(1), weights = weights)
      step_size <- 1 / iter
      beta <- beta - step_size * grad_beta
      Y_adjusted <- Y - contract_matrices(X, beta)
    }

    # Update weights
    weights <- update_weights(Y_adjusted, weights$lambda, weights$omega)
    objective_vals <- c(objective_vals, weights$val)

    # Check convergence
    if (iter >= 2) {
      improvement <- objective_vals[iter - 1] - objective_vals[iter]
      if (improvement < min_decrease^2) {
        break
      }
    }
  }

  list(
    lambda = weights$lambda,
    omega = weights$omega,
    beta = beta,
    vals = objective_vals
  )
}
