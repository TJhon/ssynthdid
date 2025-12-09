#' Setup panel data for synthetic control estimation
#'
#' @param df Input data frame
#' @param unit Unit identifier column name
#' @param time Time identifier column name
#' @param treatment Treatment indicator column name
#' @param outcome Outcome variable column name
#' @param covariates Vector of covariate column names
#' @return List with prepared data and metadata
#' @export
data_setup <- function(df, unit, time, treatment, outcome, covariates = character(0)) {

  # Convert to data.table for performance
  dt <- data.table::as.data.table(df)

  # Select and rename core columns
  core_cols <- c(unit, time, treatment, outcome)
  dt_core <- dt[, ..core_cols]
  data.table::setnames(dt_core, c("unit", "time", "treatment", "outcome"))

  # Compute treatment indicators efficiently
  dt_core[, treated := max(treatment), by = unit]
  dt_core[, ty := data.table::fifelse(treatment == 1L, time, NA_real_)]
  dt_core[, tyear := data.table::fifelse(treated == 1L, min(ty, na.rm = TRUE), NA_real_), by = unit]
  dt_core[is.na(ty), ty := 0]
  dt_core[is.na(tyear), tyear := 0]

  # Get treatment break points
  break_points <- sort(unique(dt_core$tyear))
  break_points <- break_points[break_points > 0]

  # Add covariates if present
  if (length(covariates) > 0) {
    dt_cov <- dt[, ..covariates]
    dt_core <- cbind(dt_core, dt_cov)
  }

  # Sort for efficient processing
  data.table::setorder(dt_core, treated, time, unit)

  list(
    data_ref = dt_core,
    break_points = break_points,
    covariates = covariates
  )
}


#' Convert data.table to matrix format for synthetic control
#'
#' @param dt data.table with unit, time, and value columns
#' @param value_col Name of column to extract
#' @param return_units Whether to return unit identifiers
#' @return Matrix or list with matrix and unit IDs
from_dt_to_matrix <- function(dt, value_col = "outcome", return_units = FALSE) {

  # Pivot to wide format using data.table
  dt_subset <- dt[, .(unit, time, value = get(value_col))]
  dt_wide <- data.table::dcast(dt_subset, unit ~ time, value.var = "value")

  # Extract matrix
  unit_ids <- dt_wide$unit
  mat <- as.matrix(dt_wide[, -"unit"])

  if (return_units) {
    return(list(matrix = mat, units = unit_ids))
  }

  mat
}


#' Create collapsed form of outcome matrix
#'
#' @param Y Outcome matrix (N x T)
#' @param N0 Number of control units
#' @param T0 Number of pre-treatment periods
#' @return Collapsed matrix
collapsed_form <- function(Y, N0, T0) {
  N <- nrow(Y)
  T <- ncol(Y)

  # Pre-treatment control block
  Y_11 <- Y[seq_len(N0), seq_len(T0), drop = FALSE]

  # Post-treatment control averages (column)
  Y_12 <- matrix(
    rowMeans(Y[seq_len(N0), (T0 + 1):T, drop = FALSE]),
    ncol = 1
  )

  # Pre-treatment treated averages (row)
  Y_21 <- matrix(
    colMeans(Y[(N0 + 1):N, seq_len(T0), drop = FALSE]),
    nrow = 1
  )

  # Post-treatment treated average (scalar)
  Y_22 <- mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])

  rbind(
    cbind(Y_11, Y_12),
    cbind(Y_21, Y_22)
  )
}


#' Sum normalize a vector
#'
#' @param x Numeric vector
#' @return Normalized vector that sums to 1
sum_normalize <- function(x) {
  total <- sum(x)
  if (total != 0) {
    x / total
  } else {
    rep(1 / length(x), length(x))
  }
}


#' Sparsify weight vector by setting small values to zero
#'
#' @param v Weight vector
#' @param threshold Relative threshold (default: 0.25)
#' @return Sparsified and renormalized weights
sparsify_function <- function(v, threshold = 0.25) {
  cutoff <- max(v) * threshold
  v_sparse <- ifelse(v <= cutoff, 0, v)
  v_sparse / sum(v_sparse)
}


#' Compute pairwise sum with proper handling of different lengths
#'
#' @param x First vector
#' @param y Second vector
#' @return Pairwise sum with NA handling

pairwise_sum_decreasing = function(x, y) {
  xl <- length(x)
  yl <- length(y)
  diff_vector <- xl - yl

  diff_vector_abs <- abs(diff_vector)
  diff_na <- rep(NA, diff_vector_abs)

  # ifelse(diff_vector < 0, c(a, rep(NA, diff_vector_abs)), c(b, rep(NA, diff_vector_abs)))

  if(diff_vector < 0){
    x <- append(x, diff_na)
  }else{
    y <- append(y, diff_na)
  }

  na.x = is.na(x)
  na.y = is.na(y)
  x[is.na(x)] = min(x[!na.x])
  y[is.na(y)] = min(y[!na.y])
  pairwise.sum = x + y
  pairwise.sum[na.x & na.y] = NA
  pairwise.sum
}


