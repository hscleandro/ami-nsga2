#===================================================================================================
#
# AMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm
#
# V1.0 (2019-29-10) - R Packages Dependency List
# Copyright 2019 Leandro Corrêa
#
#===================================================================================================

## Packages needed to run the AMINSGA2 tool

if("igraph" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling igraph R package...\n")
  install.packages("igraph")
}
if("scales" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling dplyr R package...\n")
  install.packages("dplyr")
}
if("magrittr" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling magrittr R package...\n")
  install.packages("magrittr")
}
if("mco" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling mco R package...\n")
  install.packages("mco")
}
if("nsga2R" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling nsga2R R package...\n")
  install.packages("nsga2R")
}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite")
}
if("rsjon" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("rsjon", quietly = TRUE))
    install.packages("rsjon")
}

## For plot network
# if("visNetwork" %in% rownames(installed.packages()) == FALSE) {
#   cat("\nInstalling R.matlab R package...\n")
#   install.packages("R.matlab")
# }