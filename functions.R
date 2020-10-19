#===================================================================================================
#
# AMINSGA2: Detect Active Modules using Multi-Objective Evolutionary Algorithm
# Publication: https://dl.acm.org/doi/10.1145/3365953.3365957
#
# V1.0 (2020-09-09)
# Author: Leandro Corrêa @leandrohsc
#
#===================================================================================================
#source("packages.R")

library(igraph)       
library(dplyr)
library(mco)
library(nsga2R)
library(igraph)
library(magrittr)
#library(visNetwork)

# --- AMINSGA-II --- #

aminsga2 <- function(G, population_size, offspring_size, number_of_generations, 
                     tournament, change_mutation_tax, deletion_mutation_tax, d_connection, 
                     variation, truehits = NULL, test = ""){
  
  cat("Start the new population...\n")
  population <- start_population(G,population_size,d_connection,variation)

  ## Set variables that will record the results ####
  mean1 <- NULL; mean2 <- NULL; st1 <- NULL; st2 <- NULL
  var1 <- NULL; var2 <- NULL; max1 <- NULL; max2 <- NULL
  number_of_springs <- NULL
  
  ## Start the evolutionary process ####
  generation <- 1; index_pop <- NULL
  
  while(generation <= number_of_generations){
    ## Evolutionary operators ####
    cat("\nCalculating Crossoover and mutation...\n")
    newpopulation <- NULL; crloop <- TRUE; eoerr <- 1; maxError <- 500

    ## choosing father and mother ####
    while(crloop){
      ## select father and mother in the population ####
      father <- V(G)$name[(which(population[select_parents(index_pop, population_size, tournament),] == 1))]
      mother <- V(G)$name[(which(population[select_parents(index_pop, population_size, tournament),] == 1))]
      ## Crossover ####
      module_child <- crossover(G, father, mother, d_connection)
      ## Mutation ####
      module_child <- mutation(G, module_child, change_mutation_tax, deletion_mutation_tax, d_connection)
      ## Identifying if the size of the module is outside of the required range
      ## attaching the new module to the new generation, if it is Ok
      if((length(module_child) > variation[1])&(length(module_child) < variation[2])){ 
        ## initialize the vector child ####
        child <- rep(0,ncol(population))
        names(child) <- colnames(population)
        ## set the child binary vector
        child[which(names(child) %in% module_child)] <- 1 
        ## join the child in the population
        newpopulation <- rbind(newpopulation,child)  
      }else{ eoerr <- eoerr + 1  }

      if(!is.null(newpopulation)){
        ## verifying if the number of children is sufficient to compose the new generation ####
        if(nrow(newpopulation) >= newgeneration){
          crloop <- FALSE
        }
      }
      if(eoerr == maxError){
        break
      }
    }
    rownames(newpopulation) <- NULL
    number_of_springs[generation] <- nrow(newpopulation)
      
    ## Combine the old population with the new population ####
    population <- rbind(population,newpopulation)
    cat("[OK]\n\n")
    
    ## Calculating the fitness functions ####
    cat("\n\nCalculating the fitness functions...\n")
    objective <- NULL
    for(i in 1:nrow(population)){
      temp <- getObjectives(population[i,],G,variation)
      objective <- rbind(objective,c(temp$za,temp$nest))
    }
    colnames(objective) <- c("za","nest")
    cat("[OK]\n")
    
    ##  Note the indicator of the objectives #### 
    mean1[generation] <- mean(objective[,1]); mean2[generation] <- mean(objective[,2])
    var1[generation] <- var(objective[,1]); var2[generation] <- var(objective[,2])
    st1[generation] <- sd(objective[,1]); st2[generation] <- sd(objective[,2])
    max1[generation] <- max(objective[,1]); max2[generation] <- max(objective[,2])
    
    ## Identifying the rank index of each Pareto Front #### 
    ## Transforming into a maximization problem
    objective <- objective * (-1)
    ranking <- fastNonDominatedSorting(objective)
    rnkIndex <- integer(nrow(population))
    i <- 1
    while (i <= length(ranking)) {
      rnkIndex[ranking[[i]]] <- i
      i <- i + 1
    }
    index_pop = matrix(data = 1:nrow(population), ncol = 1, nrow = nrow(population))
    colnames(index_pop) <- "index"
    index_pop <- cbind(index_pop,objective)
    index_pop <- cbind(index_pop,rnkIndex)
    
    ## calculating the crowding distance of each front #### 
    objRange <- apply(objective, 2, max) - apply(objective, 2, min)
    cd <- crowdingDist4frnt(index_pop,ranking,objRange)
    index_pop <- cbind(index_pop,apply(cd,1,sum))
    colnames(index_pop)[5] <- "crowdingDist"
    
    ## Choosing the best individuals per tournament #### 
    matingPool <- tournamentSelection(index_pop,nrow(population),tournament)
    matingPool <- as.data.frame(matingPool)
    matingPool <- matingPool[order(matingPool$rnkIndex),]
    matingPool <- matingPool[1:population_size,]
    
    ## Update the population 
    population <- population[matingPool$index,] 
    index_pop <- index_pop[matingPool$index,]
    index_pop[,"index"] <- 1:nrow(index_pop)
    
    ## Print informations #### 
    index <- which(matingPool$rnkIndex %in% c(1,2))
    len <- NULL;hits <- NULL
    for(i in 1:length(index)){
      temp <- length(which(population[i,] == 1))
      len <- c(len,temp)
    }
    info_table <- matingPool[index,c("za","nest","rnkIndex")]
    info_table <- cbind(1:nrow(info_table),info_table)
    info_table$za <- info_table$za * -1
    info_table$nest <- info_table$nest * -1
    rownames(info_table) <- 1:nrow(info_table)

    if(!is.null(truehits)){
      for(ix in 1:nrow(info_table)){
          res <- names(which(population[ix,] != 0))
          inter <- intersect(res,as.character(truehits))
          hits[ix] <- length(inter)
        }
      info_table <- cbind(info_table,len,hits)
    }else{
      info_table <- cbind(info_table,len)
    }
    print(head(info_table))
    cat("\n")
    cat("mean objective 1",mean1[generation],"\n")
    cat("mean objective 2",mean2[generation],"\n\n")
    
    ## Verify if the population converged to a single individual #### 
    test_pop <- as.data.frame(population)
    test_pop <- test_pop %>% distinct()
    if(nrow(test_pop) == 1){
      cat("\n--it converged!--(",generation,")\n")
      generation <- number_of_generations
    }
    
    cat("generation:",generation,"\n")
    generation <- generation + 1
  }
  
  output_list <- list()
  
  output_list[[1]] <- mean1; output_list[[2]] <- mean2; output_list[[3]] <- st1; 
  output_list[[4]] <- st2; output_list[[5]] <- max1; output_list[[6]] <- max2;
  output_list[[7]] <- var1; output_list[[8]] <- var2; output_list[[9]] <- info_table; 
  output_list[[10]] <- population; output_list[[11]] <- number_of_springs
  names(output_list) <- c("mean_obj1","mean_obj2","st_obj1","st_obj2","max_obj1","max_obj2",
                          "var_obj1","var_obj2","information","population","number_of_springs")
  return(output_list)
  
}

# --- GA Functions --- #

create_cluster <- function(ppi,d,k){
  seeds <- sample(V(ppi)$name, size = 1, replace = F)
  stopp <- 1
  while(stopp < k) {
    selectt <- names(unlist(neighborhood(ppi, order = d, nodes = seeds)))
    selectt <- setdiff(selectt,seeds)
    rnd <- sample(selectt, size = 1, replace = F)
    
    seeds <- union(seeds,rnd)
    stopp <- stopp + 1
  }
  return(seeds)
}

create_cromossom <- function(ppi,d,k = c(10,20)){
  # Just transform the clusters in a binary vector representation
  if(k[1] == k[2]){
    k1 <- k[1]
  }
  else{
    k1 <- sample(k[1]:k[2], size = 1)
  } 
  cromoss <- rep(0,length(V(ppi)))
  names(cromoss) <- V(ppi)$name
  cluster <- create_cluster(ppi,d,k1)
  index <- which(names(cromoss) %in% cluster)
  cromoss[index] <- 1
  
  return(cromoss)
}

start_population <- function(ppi,popSize=100,d,k=c(10,10)){
  # Just set a number(popSize) of cromossomes in a matrix
  population <- matrix(data = 0, ncol = length(V(ppi)$name), nrow = popSize)
  colnames(population) <- V(ppi)$name
  
  for(i in 1:popSize){
    population[i,] <- create_cromossom(ppi,d,k)
  }
  
  return(population)
}

crossover <-  function(G,P1,P2,d){
  expandSubnetworkFromSeed <- function(G,seed,P,k,d){
    subnetwork <- seed
    while(length(subnetwork) < k){
      neighbors_d <- names(unlist(neighborhood(G, order = d, nodes = subnetwork)))
      candidate <- intersect(neighbors_d,P)
      if(length(candidate) == 1){
        v <- candidate
      }else{
        v <- sample(candidate, size = 1)
      }
      subnetwork <- union(subnetwork, v)
    }
    return(subnetwork)
  }
  intersection <- intersect(P1,P2)
  if(length(intersection) == length(P1)){
    return(P1)
  }
  ## Identifying Cutoof Points and Bridges Vertices
  if(length(intersection) == 0){
    sh_path <- shortest_pathways(G,P1,P2)
    cutoff1 <- intersect(sh_path, P1)
    cutoff2 <- intersect(sh_path, P2)
    bridge <- shortest_pathways(G,cutoff1,cutoff2)
  }
  else{
    if(length(intersection) == 1){
      cutoff1 <- intersection
    }else{
      cutoff1 <- sample(intersection, size = 1)
    }
    bridge <- NULL
  }
  ## selecting the size of each part that each parent will contribute
  r <- sample(1:100, size = 1) *0.01
  size1 <- round(length(P1) * r)
  if(size1 == 0){size1 <- 1}
  size2 <- round(length(P2) * (1-r))
  if(size2 == 0){size2 <- 1}
  subnetwork1 <- expandSubnetworkFromSeed(G,cutoff1,P1,size1,d)
  if(length(intersection) > 0){
    cutoff2 <- intersect(subnetwork1,P2)
    size2 <- size2 + length(cutoff2)
    if(size2 > length(P2)){
      size2  <- length(P2) - 1
    }
  }
  
  subnetwork2 <- expandSubnetworkFromSeed(G,cutoff2,P2,size2,d)
  
  return(union(subnetwork1, union(bridge, subnetwork2)))
}

IsDconnected <- function(G,N,d){
  if(length(N) == 1){
    vqueue <- N 
  }else{
    vqueue <- sample(N, size = 1)
  }
  i <- 1
  while((i <= length(vqueue))&(length(vqueue) < length(N))){
    neighbors_d <- names(unlist(neighborhood(G, order = d, nodes = vqueue)))
    neighbors_d <- intersect(neighbors_d, N)
    for(v in neighbors_d){
      if(length(which(v %in% vqueue)) == 0){
        vqueue <- union(vqueue,v)
      }
    }
    i <- i + 1
  }
  
  if(length(vqueue) == length(N)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

mutation <- function(G, individual, p_s = 0.1, p_r = 0.08, d = 2){
  mutated <- individual
  for(v in individual){
    shrunk <- setdiff(mutated,v)
    if(IsDconnected(G,shrunk,d)){
      del <- sample(1:100, size = 1) * 0.01
      if(del <= p_r){
        mutated <- shrunk
      }else{
        flip <- sample(1:100, size = 1) * 0.01
        neighbors_d <- names(unlist(neighborhood(G, order = d, nodes = shrunk)))
        if((flip <= p_s)&(length(neighbors_d) > 0)){
          if(length(neighbors_d) == 1){
            sorted <- neighbors_d
          }
          else{
            sorted <- sample(neighbors_d, size = 1) 
          }
          mutated <- union(shrunk,sorted)
        }
      }
    }
  }
  return(mutated)
}

getObjectives <- function(cromossome, ppi, k = c(10,20)){
  ## Passing from binary format to list of genes format
  module <- V(ppi)$name[(which(cromossome == 1))]
  len <- length(module)
  ## If the module size is biggest than the user given interval the objectives return 0
  if((len < k[1])||(len > k[2])){
    nest <- 0; za <- 0
  }else{
    ## calculating the metric 'Za'-score - Ideker et al. (2002)
    index_module <- which(V(ppi)$name %in% module)
    z_score <- 1 - V(ppi)$fc[index_module] 
    
    ## Calculating inverse CDF
    z_score <- qnorm(z_score, mean = 0, sd = 1)
    za <- sum(z_score)/sqrt(len)
    
    ## calculating the mean expression of neighbors of each vertice in the module as in NEST - Jiang P., et all (2015)
    nexp <- NULL
    for(ind in 1:len){
      neighborss <- setdiff(names(unlist(neighborhood(ppi, order = 1, nodes = as.character(module[ind])))),as.character(module[ind]))
      index <- which(V(ppi)$name %in% neighborss)
      z_score <- 1 - V(ppi)$fc[index]
      
      ## Calculating inverse CDF
      z_score <- qnorm(z_score, mean = 0, sd = 1)
      nexp[ind] <- sum(z_score)/length(index)
    }
    nest <- mean(nexp)
    if(is.na(nest)){nest <- 0}
  }
  
  reply <- list()
  reply[[1]] <- za;  
  reply[[2]] <- nest; 
  reply[[3]] <- len
  names(reply) <- c("za","nest","len")
  
  return(reply)
}

select_parents <- function(index_pop, popSize, tournament){
  if(is.null(index_pop)){
    parents <- sample(1:popSize, size = 1, replace = F)
    return(parents)
  }
  ## Select by random the candidates individuals
  candidates <- sample(1:popSize, size = tournament, replace = F)
  candidates <- as.data.frame(index_pop[candidates,])
  ## Order by Front of Pareto
  ord <- order(candidates$rnkIndex)
  candidates <- candidates[ord,]
  ## Order by crowding distance
  rnk <- unique(candidates$rnkIndex)
  final_candidates <- NULL
  for(i in 1:length(rnk)){
    aux <- which(candidates$rnkIndex == rnk[i])
    index <- order(candidates$crowdingDist[aux], decreasing = T)
    aux <- candidates[aux,]
    aux <- aux[index,]
    final_candidates <- rbind(final_candidates,aux)
  }
  ## Select the best individuals and return its indexes
  parents <- final_candidates$index[1]
  
  return(parents)
}

groupingParetoSolutions <- function(result, G, d, MAX){
  info_table <- results$information
  info_table %<>%
    filter(rnkIndex == 1)
  
  set_S <- NULL
  set_Q <- NULL
  maxsize <- 0
  for(i in 1:nrow(info_table)){
    set_Q <- names(which(results$population[info_table[i,"indexes"],] == 1))
    count <- 1
    for(j in i:nrow(info_table)){
      q <- names(which(results$population[info_table[j,"indexes"],] == 1))
      if(length(union(set_Q, q)) < MAX){
        if(IsDconnected(G,union(set_Q, q),d)){
          set_Q <- union(set_Q, q)
          count <- count + 1
        }
      }
    }
    if(count > maxsize){
      set_S <- set_Q
      maxsize <- count
    }
  }
  return(set_S)
}

# --- support functions --- #

'%!in%' <- function(x,y)!('%in%'(x,y))

minMax <- function(vet, min=0, max=1){
  v_n <- NULL
  for(i in 1:length(vet)){
    v_n[i] <- min + (vet[i] - min(vet))/(max(vet) - min(vet)) *   (max-min)
    v_n[i] <- round(v_n[i], 5)
  }
  
  return(v_n)
}

shortest_pathways <- function(ppi,P1,P2){
  shortest_paths <- NULL
  for(p in P1){
    aux <- get.shortest.paths(ppi,from = as.character(p), to = as.character(P2), mode = "all")$vpath
    shortest_paths <- c(shortest_paths,aux)
  }
  min <- length(V(ppi)); index <- 0; 
  ## Calculating shortest paths between the all vertex of the two parents
  ## identifying the lowest shortest path (bridge)
  len <- NULL
  for(i in 1:length(shortest_paths)){
    len[i] <- length(names(shortest_paths[[i]]))
  }
  index <- which(len == min(len))
  index <- sample(index, size = 1)
  temp_bridge <- names(shortest_paths[[index]])
  return(temp_bridge)
}

isempty <- function(verify){
  if(length(verify) == 0){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# plot_gg <- function(gp, ind_1, ind_2){ 
#   col <- rep("#e5e5ff",length(V(gp)))
#   col[which(V(gp)$name %in% ind_1)] <- "#cf1b3c"
#   col[which(V(gp)$name %in% ind_2)] <- "#0000ff"
#   intesc <- intersect(ind_1,ind_2)
#   if(length(intesc) > 0){
#     col[which(V(gp)$name %in% intesc)] <- "#ff6600"
#   }
#   V(gp)$color <- col
#   E(gp)$color <- "silver"
# 
#   dat <- toVisNetworkData(gp)
#   nodes_path <- dat[[1]]
#   edges_path <- dat[[2]]
# 
#   visNetwork(nodes_path, edges_path,
#              main = "Interactive Network - STRING") %>%
#     visOptions(selectedBy = "color",
#                nodesIdSelection = list(enabled = TRUE), #TRUE,
#                highlightNearest = list(enabled = TRUE, degree = 1,
#                                        hover = FALSE)) %>%
#     visIgraphLayout()
# }

