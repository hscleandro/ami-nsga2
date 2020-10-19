#===================================================================================================
#
# AMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm
# Publication: https://dl.acm.org/doi/10.1145/3365953.3365957
#
# V1.0 (2020-09-09)
# Author: Leandro Corrêa @leandrohsc
#
#===================================================================================================
library(rjson)
args <- commandArgs(trailingOnly = TRUE)
LOCAL_DIR <- getwd()
if(grepl("benchmark",LOCAL_DIR)){
  LOCAL_DIR <- gsub("benchmark","\\1",LOCAL_DIR)
  TEST = TRUE
  input <- NULL
}else{ 
  LOCAL_DIR <- paste0(LOCAL_DIR,"/")
  TEST = FALSE
  if(file.exists("config.json")){
    input <- try(fromJSON(file = paste0(LOCAL_DIR,"config.json")))
    #input <- fromJSON(txt = paste0(LOCAL_DIR,"config.json"))
  }else{ input <- NULL}
}

OUTPUT_PATH <- getwd()

if("--help" %in% args){
  cat('\nAMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm\n')
  cat('\n PARAMETERS:\n')
  cat('\n --net: Edge list containg the network structure [mandatory];')
  cat('\n --profile: csv table containg in the first column the hgnc_symbol of each gene and in the second column its respectives transcript values [mandatory];')
  cat('\n --ngeneration: Parameter to determine the number of generations [integer];')
  cat('\n --popsize: Parameter to determine the size of population [integer];')
  cat('\n --minsize: Threshold parameter to determine the minimum size of the output [integer];')
  cat('\n --maxsize: Threshold parameter to determine the maximum size of the output [integer];')
  cat('\n --offspingsize Parameter to determine the size of each offspring population (based in the percent of inicial population) [integer];')
  cat('\n --changemutationtax: Parameter to determine the tax of mutation change [0..1];')
  cat('\n --delmutationtax: Parameter to determine the tax of mutation deletion [0..1];')
  cat('\n --tournament: Parameter to determine the number of candidates in each torunament selection [integer];')
  cat('\n --dconnection: Parameter to determine the number d conncetion in the aminsga2 module [integer];')
  cat('\n --output: output adress to store the results [mandatory]. \n\n')

  #cat("Warning: Don't forget to fill in the config.json file. \n\n")
  
}else{
  cat("aqui!2")
  if("--net" %in% args){
    NETWORK_PATH <- args[which(args == "--net") + 1]
  }else{
    stop("ERRO! Parameter --net is missing!")
  }
  if("--profile" %in% args){
    WEIGHT_PATH <- args[which(args == "--profile") + 1]
  }else{
    stop("ERRO! Parameter --profile is missing!")
  }
  if("--output" %in% args){
    OUTPUT_PATH <- args[which(args == "--output") + 1]
  }else{
    stop("ERRO! Parameter --output is missing!")
  }
  if("--ngeneration" %in% args){
    input$number_of_generations <- as.numeric(args[which(args == "--ngeneration") + 1])
  }else if(is.null(input$number_of_generations)){
    input$number_of_generations <- 50 #standard value
  }
  if("--popsize" %in% args){
    input$population_size <- as.numeric(args[which(args == "--popsize") + 1])
  }else if(is.null(input$population_size)){
    input$population_size <- 100 #standard value
  }
  if("--minsize" %in% args){
    input$lower_range <- as.numeric(args[which(args == "--minsize") + 1])
  }else if(is.null(input$lower_range)){
    input$lower_range <- 10 #standard value
  }
  if("--maxsize" %in% args){
    input$upper_range <- as.numeric(args[which(args == "--maxsize") + 1])
  }else if(is.null(input$upper_range)){
    input$upper_range <- 30 #standard value
  }
  if("--offspingsize" %in% args){
    input$offspring_size <- as.numeric(args[which(args == "--offspingsize") + 1])
  }else if(is.null(input$offspring_size)){
    input$offspring_size <- 0.7 #standard value
  }
  if("--changemutationtax" %in% args){
    input$change_mutation_tax <- as.numeric(args[which(args == "--changemutationtax") + 1])
  }else if(is.null(input$change_mutation_tax)){
    input$change_mutation_tax <- 0.3 #standard value
  }
  if("--delmutationtax" %in% args){
    input$deletion_mutation_tax <- as.numeric(args[which(args == "--delmutationtax") + 1])
  }else if(is.null(input$deletion_mutation_tax)){
    input$deletion_mutation_tax <- 0.08 #standard value
  }
  if("--tournament" %in% args){
    input$tournament <- as.numeric(args[which(args == "--tournament") + 1])
  }else if(is.null(input$tournament)){
    input$tournament <- 4 #standard value
  }
  if("--dconnection" %in% args){
    input$d_connection <- as.numeric(args[which(args == "--dconnection") + 1])
  }else if(is.null(input$d_connection)){
    input$d_connection <- 2 #standard value
  }
 
  source(paste0(LOCAL_DIR,"main.R"))
  
}
