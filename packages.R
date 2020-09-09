#===================================================================================================
#
# AMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm
#
# V1.0 (2019-29-10) - R Packages Dependency List
# Copyright 2019 Leandro Corrêa
#
#===================================================================================================

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
## For teste Robinson 2017 simulated dataset
# if("R.matlab" %in% rownames(installed.packages()) == FALSE) {
#   cat("\nInstalling R.matlab R package...\n")
#   install.packages("R.matlab")
# }
## For plot network
# if("visNetwork" %in% rownames(installed.packages()) == FALSE) {
#   cat("\nInstalling R.matlab R package...\n")
#   install.packages("R.matlab")
# }