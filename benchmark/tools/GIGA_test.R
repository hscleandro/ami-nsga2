#===================================================================================================
#
# GIGA test using simulated data
# This script is part of the tests used for the AMINSGA2 tool
#
# Author: Leandro Corrêa @leandrohsc
#
#===================================================================================================

#NETWORK_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_connection_list.edgelist'
#WEIGHT_PATH <- '/home/leandro/aminsga2/benchmark/output/sim2_weight.csv'
#OUTPUT_PATH <- '/home/leandro/aminsga2/benchmark/output/'

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
  cat("\n --th: Threshold parameter to determine the minimum size of the output (GIGA).\n\n")
  
  cat("Ex: Rscript GIGA_test.R -n /network.edgelist -w /weight.csv -o /output/ --th 50\n\n")
  
}else{
  ## It is important to install these libraries before run this script
  library(dplyr)
  library(magrittr)
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
    
    INPUT_PATH <- paste0(OUTPUT_PATH,"input_giga.txt")
    INPUT_NET <- paste0(OUTPUT_PATH,"graph.txt")
    TEMP_PATH <- paste0(OUTPUT_PATH,"temp.txt")
    OUTPUT_PATH <- paste0(OUTPUT_PATH,"giga_result.json")
    
  }
  if("--th" %in% args){
    arg.th <-  which(args == "--th") + 1
    TH <- as.numeric(as.character(args[arg.th]))
    #TH = 20  
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
  names(transcript) <- genesname[1:length(transcript)]
  
  ## rename PPI network
  for(i in 1:nrow(PPI)){
    PPI[i,1] <- names(transcript)[as.numeric(PPI[i,1])]
    PPI[i,2] <- names(transcript)[as.numeric(PPI[i,2])]
  }

  write.table(PPI, file = INPUT_NET,
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Obtaining the gene names
  genes <- names(transcript)
  # Ranking the genes base in the transcript level
  test.rp <- as.numeric(transcript * 10^4)
  ranks <- rank(test.rp)
  input_giga <- cbind(genes,as.numeric(test.rp),as.numeric(ranks))
  colnames(input_giga) <- c("genes", "rank.prod", "ranks")
  input_giga <- input_giga[order(ranks),]

  write.table(input_giga, file = INPUT_PATH,
            quote = FALSE, sep = "\t", row.names = FALSE)
  # Executing giga
  sys.text <- paste("perl ", GIGA_PATH, " -i", INPUT_PATH, 
                    " -n", INPUT_NET, " -Ftxt", " -t", 1, " -m", TH, " -o", TEMP_PATH, sep="");
  system(sys.text)
  
  # Extracting the results
  test <- readLines(TEMP_PATH)
  paths <- which(1:length(test)%in%grep("^-", test)==F)[1:2];
  genes <- toString(unlist(do.call(rbind,sapply(test[(paths[1] + 1):(paths[2] - 1)], FUN=function(x)strsplit(x, "\t")))[,2]))
  module <- strsplit(genes,",")[[1]]
  module <- gsub(" ", "", module, fixed = TRUE)
  
  # Converting the gene names to the names used in the simulation
  selected <- which(names(transcript) %in% module)
  selected <- weight$V1[selected]
  
  write_json(selected, path = OUTPUT_PATH)
  
}
