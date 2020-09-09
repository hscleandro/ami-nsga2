#===================================================================================================
#
# COSINE test using simulated data
# This script is part of the tests used for the AMINSGA2 tool
#
# Copyright 2020 Leandro Corrêa @leandrohsc
#
#===================================================================================================

#NETWORK_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_connection_list.edgelist'
#WEIGHT_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_weight.csv'
#OUTPUT_PATH <- '/home/leandro/aminsga2/benchmark/output/cosine_result.json'

if("dplyr" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling dplyr R package...\n")
  install.packages("dplyr")
}
if("magrittr" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling magrittr R package...\n")
  install.packages("magrittr")
}
if("igraph" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling igraph R package...\n")
  install.packages("igraph")
}
if("COSINE" %in% rownames(installed.packages()) == FALSE) {
  cat("\nInstalling COSINE R package...\n")
  install.packages("COSINE")
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
  cat("\n -o: Directory to receive the BioNet active module results;")
  cat("\n --min: The minimal size of clusters (COSINE parameter);")
  cat("\n --popSize: The number of individuals in the initial population;")
  cat("\n --int: Number of generations in the GA (COSINE parameter).\n\n")
  
  cat("Ex: Rscript COSINE_test.R -n /network.edgelist -w /weight.csv -o /output/ --min 10 --popSize 100 --int 300\n\n")
  
}else {
  ## It is important to install these libraries before run this script
  library(dplyr)
  library(magrittr)
  library(igraph)
  library(COSINE)
  library(jsonlite)
  
  minSize = 10; n_iter = 100; pop_size <- 100
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
    OUTPUT_PATH <- paste0(OUTPUT_PATH,"cosine_result.json")
  }
  if("--min" %in% args){
    arg.min <-  which(args == "--min") + 1
    minSize <- as.numeric(as.character(args[arg.min]))
  }
  if("--popSize" %in% args){
    arg.popSize <-  which(args == "--popSize") + 1
    pop_size <- as.numeric(as.character(args[arg.popSize]))
  }
  if("--int" %in% args){
    arg.int <-  which(args == "--int") + 1
    n_iter <- as.numeric(as.character(args[arg.int]))
  }
########################################################################################################

  ## upload network
  PPI <- read.csv(file = NETWORK_PATH, header = F, sep = " ")
  PPI <- as.matrix(PPI[,c(1,2)])
  # Vector python starts in 0, R starts in 1
  PPI[,1] <- PPI[,1] + 1
  PPI[,2] <- PPI[,2] + 1
  
  
  ## upload pvalues
  weight <- read.csv(file = WEIGHT_PATH, header = F)
  transcript <- weight$V2
  transcript <- 1 - as.numeric(transcript)
  names(transcript) <- 1:length(transcript)
  
  ## Fv = 1 - simulated p-values
  fv <- 1 - as.numeric(transcript)
  
  ## Calculate scaled Nodescore
  scaled_node_score <- NULL
  for(i in 1:length(transcript)){
    fv_i <- fv[i]
    ## Nodescore = Fv - mean(Fv) / stardad deviton(Fv)
    scaled_node_score[i] <- (fv_i - mean(fv)) / sd(fv)
  }
  names(scaled_node_score) <- names(transcript)
  
  ## Calculate scaled Edgescore
  edgescore <- NULL
  for(i in 1:nrow(PPI)){
    index_ey <- which(names(scaled_node_score) == as.character(PPI[i,1]))
    index_ex <- which(names(scaled_node_score) == as.character(PPI[i,2]))
    ey <- scaled_node_score[index_ey]
    ex <- scaled_node_score[index_ex]
    edgescore[i] <- (ey + ex) / 2
  }
  scaled_edge_score <- NULL
  for(i in 1:length(edgescore)){
    ecf_i <- edgescore[i]
    ## Edgescore = ECF(e) - mean(ECF) / stardad deviton(ECF)
    scaled_edge_score[i] <- (ecf_i - mean(edgescore)) / sd(edgescore)
  }
  
  ## Choice the best lambda
  klist <- minSize; 
  quantiles <- get_quantiles_PPI(scaled_node_score, scaled_edge_score, PPI, klist, pop_size)
  lambda <- quantiles[[2]]
  
  ## Exectue COSINE for the 5 best lambdas
  best_score <- 0
  for(i in 1:length(lambda)){
    lbd <- lambda[i]
    GA_result <- GA_search_PPI(lambda=lbd,scaled_node_score,scaled_edge_score,PPI,
                               num_iter=n_iter, muCh=0.05, zToR=10, minsize=minSize)
    
    bsscore <- GA_result$Best_Scores
    if(bsscore > best_score){
      ## Best COSINE results 
      subnet <- GA_result$Subnet
      best_score <- bsscore
    }
  }
  selected <- weight$V1[as.numeric(subnet)]
  
  write_json(selected, path = OUTPUT_PATH)
}