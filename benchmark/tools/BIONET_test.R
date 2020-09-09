#===================================================================================================
#
# Bionet test using simulated data
# This script is part of the tests used for the AMINSGA2 tool
#
# Copyright 2020 Leandro Corrêa @leandrohsc
#
#===================================================================================================

#NETWORK_PATH <- '/home/leandro/workspace/aminsga2/benchmark/output/sim2_connection_list.edgelist'
#WEIGHT_PATH <- '/home/leandro/workspace/aminsga2/benchmark/output/simt2_weight.csv'
#OUTPUT_PATH <- '/home/leandro/workspace/aminsga2/benchmark/output/bionet_result.json'

if("BioNet" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("BioNet")
}
if("RBGL" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("BioNet")
}
if("igraph" %in% rownames(installed.packages()) == FALSE) {
    install.packages("igraph")
}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite")
}
args <- commandArgs(trailingOnly = TRUE)

if("--help" %in% args){
  system("clear")
  cat('\nPARAMETER OPTIONS:\n')
  cat('\n -n: Edge list containg the network structure;')
  cat('\n -w: List containg weight (transcript signal) of each vertex contained in the edge list;')
  cat("\n -o: Directory to receive the BioNet active module results.")
  cat("\n --fdrs: BioNet parameter, default = 0.001.\n\n")

  cat("Ex: Rscript BIONET_test.R -n /network.edgelist -w /weight.csv -o /output/ --fdrs 0.001\n\n")
  
}else {
  ## It is important to install these libraries before run this script
  library(BioNet)
  library(jsonlite)
  library(RBGL)
  library(igraph)
  
  fdrs <- 0.001
  if("-n" %in% args){
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
    OUTPUT_PATH <- paste0(OUTPUT_PATH,"bionet_result.json")
  }
  if("--fdrs" %in% args){
    arg.fdrs <-  which(args == "--fdrs") + 1
    fdrs <- as.numeric(as.character(args[arg.fdrs]))
  }

########################################################################################################
  minMax <- function(vet, min=-1.653694, max=4.246089){
    v_n <- NULL
    for(i in 1:length(vet)){
      v_n[i] <- min + (vet[i] - min(vet))/(max(vet) - min(vet)) *   (max-min)
      v_n[i] <- round(v_n[i], 3)
    }
    return(v_n)
  }
  
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
  
  subnet <- igraph.to.graphNEL(graph.data.frame(PPI, directed=F))
  subnet <- rmSelfLoops(subnet)
  fb <- fitBumModel(transcript, plot=FALSE)
  scores <- scoreNodes(subnet, fb, fdr=fdrs)
  if(length(scores[which(scores > 0)]) == 0){
    nm <- names(scores)
    scores <- minMax(scores)
    names(scores) <- nm 
  }
  module <- try(runFastHeinz(subnet, scores))
  if("try-error" %in% class(module)){
    selected <- NULL
  } else{
    selected <- module@nodes 
  }
  selected <- weight$V1[as.numeric(selected)]
  
  write_json(selected, path = OUTPUT_PATH)
} 
