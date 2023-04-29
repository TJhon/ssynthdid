#' California proposition 99
#'
#' A dataset containing per-capita cigarette consumption (in packs).
#' In year 1989 California imposed a Tobacco tax. The column `treated` is 1 from then on for California.
#'
#' @format A data frame with 1209 rows and 4 variables:
#' \describe{
#'   \item{State}{US state name, character string}
#'   \item{Year}{Year, integer}
#'   \item{PacksPerCapita}{per-capita cigarette consumption, numeric}
#'   \item{treated}{the treatmed indicator 0: control, 1: treated, numeric}
#' }
#' @return california_prop99
#' @export
california <- function(){
  readr::read_delim(
    "https://github.com/d2cml-ai/Synthdid.jl/raw/stag_treat/data/california_prop99.csv",
    show_col_types = F)
}


#' Quota Data set
#'
#' @export
quota <- function(){
  readr::read_csv("https://github.com/d2cml-ai/Synthdid.jl/raw/stag_treat/data/quota.csv", show_col_types = F)
}
#
