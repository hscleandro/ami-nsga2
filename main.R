#===================================================================================================
#
# AMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm
# Publication: https://dl.acm.org/doi/10.1145/3365953.3365957
#
# V1.0 (2020-09-09) 
# Author: Leandro Corrêa @leandrohsc
#
#===================================================================================================
library(jsonlite)

#LOCAL_DIR = '~/aminsga2/'
#input <- try(fromJSON(file = paste0(LOCAL_DIR,"config.json")))

#OUTPUT_PATH <- '~/aminsga2/benchmark/output/'
#NETWORK_PATH <- '~/aminsga2/benchmark/output/sim1_connection_list.edgelist'
#WEIGHT_PATH <- '~/aminsga2/benchmark/output/sim1_weight.csv'
#HIT_PATH <- '~/aminsga2/benchmark/output/sim1_hittrue.csv'

#TEST = TRUE

source(paste0(LOCAL_DIR,"functions.R"))

if(input$lower_range > input$upper_range){
  cat("\n\n")
  cat("STOP! minLimit can not be smaller than maxLimit! \n")
  stop("Error", call. = FALSE)
}else{
  variation = NULL
  variation[1] = input$lower_range
  variation[2] = input$upper_range
}

## upload network
PPI <- read.csv(file = NETWORK_PATH, header = F, sep = " ")
PPI <- as.matrix(PPI[,c(1,2)])

## upload p-values
weight <- read.csv(file = WEIGHT_PATH, header = F)
profile <- weight
colnames(profile) <- c("names","values")

if(TEST){
  cat("\n Benchmark test.\n")
  # Vector python starts in 0, R starts in 1
  PPI[,1] <- PPI[,1] + 1
  PPI[,2] <- PPI[,2] + 1
  profile$names <- 1:length(weight$V1)
  profile$values <- 1 - profile$values 
}

## Normalizing the values of the expression profile (if they are outside the 0..1 range)
if((max(profile$values) > 1)||(min(profile$values) < 0)){
  profile$values <- minMax(profile$values, min = 1e-21, max = 0.999999)
}

# starting PPI igraph network
ppi <- graph.edgelist(PPI)
V(ppi)$name <- profile$names

## Using high score values to initialize the weight of the network vertices.
maxx <- max(profile$values)
## Indexing the pvalue from transcript analysis 
V(ppi)$fc <- rep(maxx,length(V(ppi)$name))
for(i in 1:nrow(PPI)){
  index <- which(V(ppi)$name[i] == as.character(profile$names))
  if(length(index) > 0 ){
    V(ppi)$fc[i] <- as.numeric(profile$values[index])
  }
}

cat("\n\n## Executing AMINSGA-II ## \n\n")

cat("Number of population:", input$population_size,"\n")
newgeneration <- trunc(input$population_size * input$offspring_size)
cat("Number of individuals per generation:",newgeneration,"\n")
cat("Number of generations: ",input$number_of_generations,"\n")
cat("Tournament: ",input$tournament,"\n")
cat("Mutation change rate: ",input$change_mutation_tax,"\n")
cat("Mutation deletion rate: ",input$deletion_mutation_tax,"\n")
cat("Minimum active module size: ",input$lower_range,"\n")
cat("Maximum active module size: ",input$upper_range,"\n")
cat("d connextion size: ",input$d_connection,"\n")
cat("[OK]\n\n")

## Calling AMINSGA function
results <- aminsga2(G = ppi, 
                    population_size = input$population_size, 
                    offspring_size = input$offspring_size, 
                    number_of_generations = input$number_of_generations, 
                    tournament = input$tournament, 
                    change_mutation_tax = input$change_mutation_tax,
                    deletion_mutation_tax = input$deletion_mutation_tax, 
                    d_connection = input$d_connection, 
                    variation = variation)

cat("\n\n## Identifying the best grouping Pareto combination.\n")
final_module <- groupingParetoSolutions(results, ppi, input$d_connection, variation[2])

if(TEST){
  selected <- weight$V1[as.numeric(final_module)]
}else{
  selected <- final_module
}

write_json(selected, path = paste0(OUTPUT_PATH,"aminsga2_result.json"))

#hits <- read.csv(HIT_PATH, header = F)
#true_positives <- intersect(as.character(hits$V1),final_module)
