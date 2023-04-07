librarian::shelf(here, purrr)

ref_code <- c("get_data.r", "solver.r", "utils.r")

here("R", ref_code) |> map(source)



staggered_synthestimate <- function(
    data, unit, time, treated, outcome, covariates = c(),  weights_sdid = list(lambda=NULL, omega=NULL),
    update_lambda = is.null(weights_sdid$lambda), update_omega = is.null(weights_sdid$omega),
    cov_method = "optimized", lambda_intercept = T, omega_intercept = T, max_iter_pre_sparsify = 100,
    max_iter = 1e4, sparsify = sparsify_function
                                    ){
  data_ref <- data_setup(data, unit, time, treated, outcome, covariates)
  tdf <- data_ref$data_ref
  break_points <- data_ref$break_points

  synthd_att_search <- function(tdf, time_eval, covariates = data_ref$covariates, cov_method = "optimized"){
    df_y <- tdf |> dplyr::filter(tyear %in% c(0, time_eval))

    N1 <- df_y |> dplyr::filter(treated == 1) |> dplyr::pull(unit) |> unique() |> length()
    T1 <- (tdf |> dplyr::pull(time) |> max()) - time_eval + 1
    tau_hat_wt <-  N1 * T1

    from_tdf_matrix <- function(values_wider = "outcome", df = df_y) {
      vars <- c("time", "unit", values_wider)
      df_y |> dplyr::select(dplyr::all_of(vars)) |> tidyr::pivot_wider(names_from = time, values_from = values_wider) |>
        dplyr::select(!unit) |> as.matrix()
    }

    # Y = df_y |> dplyr::select(time, unit, outcome) |> tidyr::pivot_wider(names_from = time, values_from = outcome) |>
      # dplyr::select(!unit) |> as.matrix()
    Y <- from_tdf_matrix()

    N_y <- nrow(Y)
    T_y <- ncol(Y)
    N0 <- N_y - N1
    T0 <- T_y - T1

    # Params from synthdid_estimate

    noise_level <- sd(apply(Y[1:N0, 1:T0], 1, diff))
    eta_omega <- ((N_y - N0) * (T_y - T0))^(1/4)
    eta_lambda <-  1e-6
    zeta_omega <- eta_omega * noise_level
    zeta_lambda <- eta_lambda * noise_level
    min_decrease <- 1e-5 * noise_level

    Yc <- collapsed_form(Y, N0, T0)

    ### No covariates
    if(is.null(covariates) || cov_method == "projected"){
      weights_sdid <- list()
      weights_sdid$vals <- weights_sdid$lambda_vals <- weights_sdid$omega_vals <- NULL

      Al <- Yc[1:N0, ]
      Ao <- t(Yc[, 1:T0])
      if(update_lambda){
        lambda_opt <- sc.weight.fw(
          Al, zeta = zeta_lambda, intercept = lambda_intercept, lambda=NULL,
          min.decrease = min_decrease, max.iter = max_iter_pre_sparsify)
        if(!is.null(sparsify)){
          lambda_opt <- sc.weight.fw(
            Al, zeta = zeta_lambda, intercept = lambda_intercept, lambda=sparsify(lambda_opt$lambda),
            min.decrease = min_decrease, max.iter = max_iter)
        }
        weights_sdid$lambda <- lambda_opt$lambda
        weights_sdid$lambda_vals <- weights_sdid$vals <- lambda_opt$vals
      }
      if(update_omega){
        omega_opt <- sc.weight.fw(
          Ao, zeta = zeta_omega, intercept = omega_intercept, lambda=NULL,
          min.decrease = min_decrease, max.iter = max_iter_pre_sparsify)
        if(!is.null(sparsify)){
          omega_opt <- sc.weight.fw(
            Ao, zeta = zeta_omega, intercept = omega_intercept, lambda=sparsify(omega_opt$lambda),
            min.decrease = min_decrease, max.iter = max_iter)
        }
        weights_sdid$omega <- omega_opt$lambda
        weights_sdid$omega_vals <- omega_opt$vals
        if(is.null(weights_sdid$vals)) {weights_sdid$vals <- omega_opt$vals}
        else{weights_sdid$vals <- pairwise_sum_decreasing(weights_sdid$vals, omega_opt$vals)}
      }
      X_beta <- 0
    } else if(length(covariates) > 0)
      {
      X <- map(covariates, from_tdf_matrix)
      Xc <- X |> map(collapsed_form, N0, T0)
      # print(X_covariates |> map(dim) )
      # print(dim(Yc))
      weights_sdid = sc.weight.fw.covariates(
        Yc, Xc, zeta.lambda = zeta_lambda, zeta.omega = zeta_omega,
        lambda.intercept = lambda_intercept, omega.intercept = omega_intercept,
        min.decrease = min_decrease, max.iter = max_iter,
        lambda = weights_sdid$lambda, omega = weights_sdid$omega, update.lambda = update_lambda, update.omega = update_omega)
      info_beta <- weights_sdid$beta
      names(info_beta) <- paste0("beta_", covariates)
      info_beta <- as.matrix(info_beta) |> t() |> tibble::as_tibble(.name_repair = "unique")
      X_beta <- contract3(X, weights_sdid$beta)
      # print(dplyr::glimpse(weights_sdid))
    }

    tau <- t(c(-weights_sdid$omega, rep(1 / N1, N1))) %*% (Y - X_beta) %*% c(-weights_sdid$lambda, rep(1 / T1, T1))

    info_att <- tibble::tibble(
      "time" = time_eval,
      "tau" = as.double(tau), "tau_wt" = tau_hat_wt, "N0" = N0, "T0" = T0, "N1" = N1, "T1" = T1,
      "weights_sdid" = list(weights_sdid)
      )


    if(!is.null(covariates)){
      info_att <- info_att |> dplyr::bind_cols(info_beta)
    }
    return(info_att)
  }

  att_table <- purrr::map_df(break_points, synthd_att_search, tdf = tdf)
  att_table <- att_table |> dplyr::mutate(
    weighted_tau = tau * tau_wt / sum(tau_wt)
  ) |> dplyr::relocate(weighted_tau, .after = tau_wt)
  att <- sum(att_table$weighted_tau)
  return(list(att_estimate = att, att_table = att_table))
}

staggered_synthestimate(quota() |> tidyr::drop_na(lngdp), "country", "year", "quota", "womparl", "lngdp")


staggered_synthestimate(california(), "State", "Year", "treated", "PacksPerCapita")



