#' Plots treated and synthetic control trajectories and overlays a 2x2 diff-in-diff diagram of our estimator.
#' Requires ggplot2
#' @param estimate A Ssynthdid model
#'
#' @return A list with all the plots for the different treatment periods.
#' @export
ssynthdid_plot <- function(estimate){
  T1s <- estimate$att_table$T1
  T0s <- estimate$att_table$T0
  N1s <- estimate$att_table$N1
  N0s <- estimate$att_table$N0

  omega_wgts <- estimate$att_table$weights_sdid |> purrr::map(purrr::pluck, "omega")
  lambda_wgts <- estimate$att_table$weights_sdid |> purrr::map(purrr::pluck, "lambda")

  Y_setups <- estimate$att_table$Y_beta

  times <- estimate$att_table$time

  id_time <- 1

  make_plot <- function(id_time){
    omega_hat <- omega_wgts[[id_time]]
    lambda_hat <- lambda_wgts[[id_time]]
    N0 <- N0s[id_time]
    T0 <- T0s[id_time]
    N1 <- N1s[id_time]
    T1 <- T1s[id_time]

    Y_time <- Y_setups[[1]]
    Y_c <- Y_time[1:N0, ]

    time_span <- colnames(Y_time) |> as.numeric()

    # treatment
    Y_sdid_traj <- as.vector(omega_hat %*% Y_c)
    # contorl
    Y_t <- Y_time[N0:nrow(Y_time), ] |> colMeans() |> as.vector()


    traj_value <- c(Y_t, Y_sdid_traj)

    y_min <- min(traj_value)
    y_max <- max(traj_value)
    plot_h <- y_max - y_min
    base_plot <- y_min - plot_h / 5

    colors_st <- c("#389826", '#cb3c33')

    ribbond <- tibble::tibble(t_time = time_span[1:T0], line = lambda_hat * plot_h / 3 + base_plot)

    lines <- tibble::tibble(
      t_time = time_span, Control = Y_t,
      Treatment = Y_sdid_traj
      ) |>
      tidyr::pivot_longer(!t_time)

    time <- times[id_time]


    p <-
      ggplot(lines) +
      aes(t_time) +
      geom_line(aes(y = value, color = name), size = 1) +
      labs(x = "Time", y = "", color = "", title = paste0("Time adoption: ", time)) +
      geom_vline(xintercept = time, linetype = "dashed") +
      geom_ribbon(aes(ymax = line), ymin = base_plot, data = ribbond, fill = "#9558b2") +
      scale_color_manual(
        values = colors_st, labels = c("Control", "Treatment")
      ) +
      theme_bw()

    return(p)
  }

  plots_sdid <- purrr::map(1:length(times), make_plot, .progress = T)
  names(plots_sdid) <- paste0("time_", times)
  return(plots_sdid)
}

#' Plots unit by unit difference-in-differences. Dot size indicates the weights omega_i
#' used in the average that yields our treatment effect estimate.
#'
#' @param estimate A Ssynthdid model
#'
#' @return A list with all the plots for the different treatment periods.
#' @export
ssynthdid_units_plot <- function(estimate){


  N0s <- estimate$att_table$N0; N1s <- estimate$att_table$N1
  T0s <- estimate$att_table$T0; T1s <- estimate$att_table$T1

  wgts <- estimate$att_table$weights_sdid
  omega_wgts <- wgts |> purrr::map(purrr::pluck, "omega")
  lambda_wgts <- wgts |> purrr::map(purrr::pluck, "lambda")
  Y_units <- estimate$att_table$Units
  Y_setups <- estimate$att_table$Y_beta
  att <- estimate$att_estimate |> round(2)
  atts <- estimate$att_table$tau |> round(2)
  times <- estimate$att_table$time

  make_plot <- function(id_time){
    lambda_hat <- lambda_wgts[[id_time]]; omega_hat <- omega_wgts[[id_time]]
    T1 = T1s[id_time]; T0 = T0s[id_time]; N0 = N0s[id_time]; N1 = N1s[id_time]
    yunits <- Y_units[[id_time]]
    lambda_pre <- c(lambda_hat, rep(0, T1))
    lambda_post <- c(rep(0, T0), rep(1 / T1, T1))
    omega_control <- c(omega_hat, rep(0, N1))
    omega_treat <- c(rep(0, N0), rep(1 / N1, N1))
    Y <- Y_setups[[id_time]]
    att_tau <- atts[id_time]
    text_subtitle <- ifelse(length(atts) < 2, paste0("ATT Stagered: ", att),paste0("ATT Stagered: ", att, paste("\nATT", times[id_time], ":", atts[id_time])))


    difs = as.vector(t(omega_treat) %*% Y %*% (lambda_post - lambda_pre)) - as.vector(Y[1:N0, ] %*% (lambda_post - lambda_pre))

    size_dot <- omega_hat / max(omega_hat) * 10

    weight_df <-
      tibble::tibble(
        dot_value = difs, size_dot = size_dot, units = yunits[1:N0]
        ) |>
      dplyr::mutate(
        shape_dot = ifelse(size_dot == 0, "0", "1") ,
        units = forcats::fct_reorder(units, shape_dot)
        )

    plot <- weight_df |>
      ggplot() +
      aes(x = units, y= dot_value, size = size_dot, shape = shape_dot) +
      geom_hline(yintercept = att, color = "#0b4d81", size = 1, linetype = "dotdash") +
      geom_point() +
      scale_shape_manual(values = c(4, 16)) +
      labs(
        x = "", y = "", title=paste0("Time adtoption: ", times[id_time])
        , subtitle = text_subtitle
           ) +
      {if(max(difs) > 0 & min(difs) < 0){
        geom_hline(yintercept = 0, color = "red", linetype = "dashed")
      }} +
      {if(length(atts) > 2){
       geom_hline(yintercept = atts[id_time], color = "#0b4d81", linetype = "dotted")
      }
      } +
      theme_light() +
      theme(
        legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    return(plot)
  }
  plots <- purrr::map(1:length(times), make_plot)
  names(plots) <- paste0("time_", times)
  return(plots)
}

