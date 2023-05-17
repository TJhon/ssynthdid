# dir("R", full.names = T) |> purrr::map(source)
library(Ssynthdid)
quota_df <- quota() |> tidyr::drop_na(lngdp)

a <- ssynth_estimate(quota_df, "country", "year", "quota", "womparl", covariates = c("lngdp"))
# ssynth_estimate(quota_df, "country", "year", "quota", "womparl")


data <- quota() |> tidyr::drop_na(lngdp)
unit <- "country"
time <- "year"
treated <- "quota"
outcome <- "womparl"

covariates <- c("lngdp")

weights_sdid = list(lambda=NULL, omega=NULL)
update_lambda = is.null(weights_sdid$lambda)
update_omega = is.null(weights_sdid$omega)
cov_method = "optimized"
lambda_intercept = T
omega_intercept = T
max_iter_pre_sparsify = 100
max_iter = 1e4
sparsify = sparsify_function





data_ref <- data_setup(data, unit, time, treated, outcome, covariates)
tdf <- data_ref$data_ref
break_points <- data_ref$break_points

# time_eval <- breackpoitn

time_eval <- break_points[1]


stt <- function(tdf, time_eval, covariates = data_ref$covariates){


df_y <- tdf |> dplyr::filter(tyear %in% c(0, time_eval))

N1 <- df_y |> dplyr::filter(treated == 1) |> dplyr::pull(unit) |> unique() |> length()
T1 <- (tdf |> dplyr::pull(time) |> max()) - time_eval + 1
tau_hat_wt <-  N1 * T1


from_tdf_matrix <- function(values_wider = "outcome", df = df_y, units = F) {
  vars <- c("time", "unit", values_wider)
  dataframe <-
    df_y |> dplyr::select(dplyr::all_of(vars)) |> tidyr::pivot_wider(names_from = time, values_from = values_wider)
  matrix <- dataframe |>
    dplyr::select(!unit) |> as.matrix()
  uniqid <- dataframe |> dplyr::pull(unit) |> unique()
  if(units){
    return(list(matrix, uniqid))
  }
  return(matrix)
}


y_omega <- from_tdf_matrix(units = T)
Y <- y_omega[[1]]
units_y <- y_omega[[2]]

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




X <- purrr::map(covariates, from_tdf_matrix)
Xc <- X |> purrr::map(collapsed_form, N0, T0)
weights_sdid = sc.weight.fw.covariates(
  Yc, Xc, zeta.lambda = zeta_lambda, zeta.omega = zeta_omega,
  lambda.intercept = lambda_intercept, omega.intercept = omega_intercept,
  min.decrease = min_decrease, max.iter = max_iter,
  lambda = weights_sdid$lambda, omega = weights_sdid$omega,
  update.lambda = update_lambda, update.omega = update_omega)


info_beta <- weights_sdid$beta
names(info_beta) <- paste0("beta_", covariates)
info_beta <- as.matrix(info_beta) |> t() |> tibble::as_tibble(.name_repair = "unique")
X_beta <- contract3(X, weights_sdid$beta)

Y_beta <- Y - X_beta
tau <- t(c(-weights_sdid$omega, rep(1 / N1, N1))) %*% (Y_beta) %*% c(-weights_sdid$lambda, rep(1 / T1, T1))

info_att <- tibble::tibble(
  "time" = time_eval,
  "tau" = as.double(tau), "tau_wt" = tau_hat_wt, "N0" = N0, "T0" = T0, "N1" = N1, "T1" = T1,
  "weights_sdid" = list(weights_sdid),
  "Y_beta" = list(Y_beta),
  "Units" = list(units_y)
)


if(!is.null(covariates)){
  info_att <- info_att |> dplyr::bind_cols(info_beta)
}
return(info_att)
}


stt(tdf, break_points[3])

purrr::map_df(break_points, stt, tdf = tdf, .progress = T)

ssynth_estimate(quota_df, unit, time, treated, outcome, covariates)
