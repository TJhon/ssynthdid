collapsed_form = function(Y, N0, T0) {
  N = nrow(Y); T = ncol(Y)
  rbind(cbind(Y[1:N0, 1:T0, drop = FALSE], rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])),
        cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])))
}
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

sparsify_function <- function(v){
  v = ifelse(v <= max(v) / 4, 0, v)
  return(v / sum(v))
}

# panel matrices
data_setup <- function(
    df, unit, time, treatment, outcome, covariates=c()
    ){

    data_ref <- df |>
      dplyr::select(dplyr::all_of(c(unit, time, treatment, outcome))) |>
      dplyr::rename(unit = 1, time = 2, treatment = 3, outcome = 4) |>
      dplyr::group_by(unit) |>
      dplyr::mutate(
        treated = max(treatment),
        ty = ifelse(treatment == 1, time, NA),
        tyear = ifelse(treated == 1, min(ty, na.rm = T), NA)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        dplyr::across(c(ty, tyear), tidyr::replace_na, 0)
      )
    break_points <- unique(data_ref$tyear) |> sort()

    if(length(covariates) > 0){
      df_cov <- df |> dplyr::select(dplyr::all_of(covariates))
      data_ref <- data_ref |> dplyr::bind_cols(df_cov)
    }

    data_ref <- data_ref |> dplyr::arrange(treated, time, unit)

    info_data <- list(
      "data_ref" = data_ref, "break_points" = break_points[-1], covariates = covariates
    )

    return(info_data)
}




