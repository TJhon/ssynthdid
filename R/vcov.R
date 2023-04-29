

jackknife_se <- function(estimate) {
  data_ref <- estimate$data_ref
  time_break <- estimate$att_table$time
  uniqID <- unique(data_ref$unit)

  if(any(estimate$att_table$N1 < 1)){
    stop("Jackknife standard error needs at least two treated units for each treatment period")
  }

  N <- length(uniqID)

  total_iter <- expand.grid(id = 1:N, 1:length(time_break))

  weights <- estimate$att_table$weights_sdid
  lambda_aux <- weights |> purrr::map(purrr::pluck, "lambda")
  omega_aux <- weights |> purrr::map(purrr::pluck, "omega")


  individuos <- total_iter$id
  times <- total_iter$Var2

  # ind: units, id, time_breaks

  theta_jk <- function(ind, id){
    omega_aux <- (omega_aux |> purrr::pluck(id))[-ind] |> sum_normalize()
    lambda_aux <- lambda_aux |> purrr::pluck(id)
    data_aux <- data_ref |>
      dplyr::filter(unit != uniqID[ind]) |>
      dplyr::filter(tyear %in% c(0, time_break[id]))

    Y_aux <-
      data_aux |>
      dplyr::select(time, unit, outcome) |>
      tidyr::pivot_wider(names_from = time, values_from = outcome) |>
      dplyr::select(!unit) |> as.matrix()

    yng <-  data_aux |> dplyr::count(unit) |> nrow()
    ynt <- data_aux |> dplyr::count(unit) |> dplyr::pull(n) |> min()
    N1 <- sum(data_aux$tyear == time_break[id]) / ynt

    npre <- data_aux |> dplyr::filter(time < time_break[id]) |> dplyr::count(time) |> nrow()
    T1 <- ynt - npre
    tau_wt_aux <- N1 * T1
    tau_aux <- att_mult(Y_aux, omega_aux, lambda_aux, N1, T1)

    tibble::tibble(unit = uniqID[ind], time = time_break[id], tau_aux = tau_aux[1], tau_wt_aux)
  }


  att_table <- purrr::map2_df(individuos, times, theta_jk, .progress = T)

  result_att <-
    att_table |> dplyr::arrange(unit) |> dplyr::distinct() |> dplyr::group_by(unit) |>
    dplyr::mutate(tau_wt = tau_wt_aux / sum(tau_wt_aux), att_aux = tau_aux * tau_wt) |>
    dplyr::summarise(att_aux = sum(att_aux)) |>
    dplyr::ungroup()

  se_jackknife <- (((N-1)/N) * (N - 1) * var(result_att$att_aux)) |> sqrt()
  return(list(result_att, se_jackknife, att_table))
}

bootstrap_se <- function(estimate, n_reps = 50){
  data_ref <- estimate$data_ref
  uniqueID <- data_ref |> dplyr::pull(unit) |> unique()

  theta_bt <- function(n){
    sample_unit <- sample(uniqueID, replace = T)

    eval_df <- data_ref |> dplyr::filter(unit %in% sample_unit)
    if(length(unique(eval_df$treatment)) != 2){
      theta_bt()
    }

    sample_to_df <- function(id){
      data_ref |> dplyr::filter(unit == sample_unit[id]) |> dplyr::mutate(unit = paste0(unit, id))
    }

    sampled_df <- purrr::map_df(1:length(sample_unit), sample_to_df)

    ssynth_estimate(sampled_df, "unit", "time", "treatment", "outcome")$att_estimate
  }

  att_bt = purrr::map_dbl(1:n_reps, theta_bt, .progress = T)

  se_bootstrap <- (1 / n_reps * sum((att_bt - sum(att_bt / n_reps, na.rm = T)) ^ 2, na.rm = T)) |> sqrt()
  return(se_bootstrap)
}

placebo_se <- function(estimate, n_reps = 50){
  data_ref <- estimate$data_ref
  tr_years <-
    data_ref |>
    dplyr::filter(
      time == tyear, tyear != 0
    ) |> dplyr::pull(time)
  N_tr <- length(tr_years)

  df_co <- data_ref |>
    dplyr::filter(treated == 0)
  N_co <- unique(df_co$unit) |> length()
  N_aux <- N_co -  N_tr

  theta_pb <- function(i){
    placebo_tibble <- tibble::tibble(unit = sample(unique(df_co$unit), N_tr), tyear1 = tr_years)

    aux_data <-
      df_co |>
      dplyr::full_join(placebo_tibble, by = "unit") |>
      dplyr::mutate(tyear = ifelse(is.na(tyear1), tyear, tyear1))

    aux_data <- aux_data |>
      dplyr::mutate(
        treatment = ifelse(
          tyear != 0 & time == tyear, 1, 0
        )
      ) |>
      # count(treatment)
      dplyr::group_by(unit) |>
      dplyr::mutate(tunit = max(treatment)) |>
      dplyr::ungroup()

    att <- ssynth_estimate(aux_data, "unit", "time", "treatment", "outcome")$att_estimate
    return(att)
  }


  att_pb <- purrr::map_dbl(1:n_reps, theta_pb, .progress = T)

  placebo_se <- (1 / n_reps * sum((att_pb - sum(att_pb / n_reps, na.rm = T)) ^ 2, na.rm = T)) |> sqrt()
  return(placebo_se)
}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates in the context of staggered adoption based on the following three options.
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#' @param estimate A Ssynthdid model
#' @param method The CI method. The default is "placebo"
#' @param n_reps the number of replicacions c("placebo", "bootstrap")
#'
#' @export
ss_vcov <- function(estimate, method = "placebo", n_reps){
  if(method == "placebo"){
    return(placebo_se(estimate, n_reps))
  }else if (method == "bootstrap"){
    return(bootstrap_se(estimate, n_reps))
  }else{
    return(jackknife_se(estimate))
  }
}
