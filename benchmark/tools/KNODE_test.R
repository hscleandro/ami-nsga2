#===================================================================================================
#
# KNODE test using simulated data
# This script is part of the tests used for the AMINSGA2 tool
#
# Author: Leandro Corrêa @leandrohsc
#
#===================================================================================================
if("dplyr" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling dplyr R package...\n")
  install.packages("dplyr")
}
if("SANTA" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(BiocManager)
  cat("\nInstalling SANTA R package...\n")
  BiocManager::install("SANTA")
}
if("igraph" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling igraph R package...\n")
  install.packages("igraph")
}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite")
}

#NETWORK_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_connection_list.edgelist'
#WEIGHT_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_weight.csv'
#OUTPUT_PATH <- '/home/leandro/aminsga2/benchmark/output/'
#minsize = 20

args <- commandArgs(trailingOnly = TRUE)

if("--help" %in% args){
  system("clear")
  cat('\nPARAMETER OPTIONS:\n')
  cat('\n -n: Edge list containg the network structure;')
  cat('\n -w: List containg weight (transcript signal) of each vertex contained in the edge list;')
  cat("\n -o: Directory to receive the BioNet active module results;")
  cat("\n --minsize: Threshold parameter to determine the minimum size of the output.\n\n")
  
  cat("Ex: Rscript GIGA_test.R -n /network.edgelist -w /weight.csv -o /output/ --th 50\n\n")
  
}else{
  ## It is important to install these libraries before run this script
  library(dplyr)
  library(SANTA)
  library(igraph) 
  library(jsonlite)
  
  GIGA_PATH <- paste0(getwd(),"/tools/GiGA.pl")
  GENE_PATH <- paste0(getwd(),"/tools/genes_hugo.csv")
  
  genesname <- read.csv(file = GENE_PATH)
  genesname <- as.character(genesname$V2)
  
  if("-n" %in% args){
    arg.files <-  which(args == "-n") + 1
    arg.net <-  which(args == "-n") + 1
    NETWORK_PATH <- as.character(args[arg.net])
  }
  if("-w" %in% args){
    arg.weight <-  which(args == "-w") + 1
    WEIGHT_PATH <- as.character(args[arg.weight])
  }
  if("-o" %in% args){
    arg.output <-  which(args == "-o") + 1
    OUTPUT_PATH <- as.character(args[arg.output])
    OUTPUT_PATH <- paste0(OUTPUT_PATH,"knode_result.json")
  }
  if("--th" %in% args){
    arg.th <-  which(args == "--th") + 1
    minsize <- as.numeric(as.character(args[arg.th]))
  }
  
  ########################################################################################################
  ## upload network
  PPI <- read.csv(file = NETWORK_PATH, header = F, sep = " ")
  PPI <- as.matrix(PPI[,c(1,2)])
  # Vector python starts in 0, R starts in 1
  PPI[,1] <- PPI[,1] + 1
  PPI[,2] <- PPI[,2] + 1
  
  ## upload p-values
  weight <- read.csv(file = WEIGHT_PATH, header = F)
  transcript <- weight$V2
  transcript <- 1 - as.numeric(transcript)
  names(transcript) <- 1:length(transcript)

  ppi <- graph.edgelist(PPI)
  ppi <- set.vertex.attribute(ppi, name="pheno", value=-log10(transcript)) 
  knode.results <- Knode(ppi, dist.method="diffusion", 
                         vertex.attr="pheno", verbose=FALSE)
  
  selected <- weight$V1[as.numeric(names(knode.results[1:minsize]))]
  write_json(selected, path = OUTPUT_PATH)
}
