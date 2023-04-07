library(readr)

california <- function(){
  read_delim(
    "https://github.com/d2cml-ai/Synthdid.jl/raw/stag_treat/data/california_prop99.csv",
    show_col_types = F)
}

quota <- function(){
  read_csv("https://github.com/d2cml-ai/Synthdid.jl/raw/stag_treat/data/quota.csv", show_col_types = F)
}
#
