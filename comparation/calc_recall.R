#===================================================================================================
#
# AMI-NSGAII: New evolutionary operators for active modules identification
#
# AMI-NSGAII v0.1 (2019-07-22)  
# Copyright 2019 Leandro Corrêa
#
# This file is part of AMI-NSGAII.
#
# This script compares the recall of the results to the tools:
# - AMI-NSGAII
# - Knode  (based in Robinson (2017) analisys)
# - Bionet (based in Robinson (2017) analisys)
# - COSINE
#
# The simulated data used are based on the work of Robinson (2017)
# You can get the data in .mat format at this address: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870666/bin/btx244_supp.zip
#
# REFERENCES
#
# Sean Robinson, Jaakko Nevalainen, Guillaume Pinna, Anna Campalans, J. 
# PabloRadicella, and Laurent Guyon. 2017.  Incorporating interaction 
# networks intothe determination of functionally related hit genes in 
# genomic experimentswith Markov random fields.Bioinformatics (Oxford, 
# England)33, 14 (jul 2017),i170–i179
#===================================================================================================

library(R.matlab)     # to upload the Robinson(2017) dataset
library(dplyr)
library(magrittr)
library(igraph)  
library(COSINE) 

## Address of simulated Robinson data in .mat
SIMULATED_DATA <- '/home/leandro/Data/Robinson/btx244-suppl_data/Robinson.121.sup.2/supplementary_material_code_and_data/simulateddata.mat'
## Directory indicating the address of the modules obtained by the AMI-NSGAII tool
FILES_PATH <- '/home/leandro/Data/nsga2_module3/table_module/'
## output file of the recall results of each tool
OUTPUT_PATH <- '/home/leandro/Data/nsga2_module3/COSINE/final_results_test3.csv'

data_sim <- readMat(SIMULATED_DATA)
files <- list.files(FILES_PATH)

## Ordering the modules files
nmod <- NULL
for(i in 1:length(files)){
  aux <- strsplit(as.character(files[i]),"_")[[1]][1]
  aux <- gsub("[net]","\\1",aux)
  nmod[i] <- as.numeric(aux)
}
index <- order(nmod, decreasing = F)
files <- files[index]

## Getting the best ami-nsgaII modules
k <- 1; size <- NULL; final_table <- NULL
for(f in files){
  NEW_PATH <- paste0(FILES_PATH,f)
  front_pareto <- read.csv(NEW_PATH)
  ## Identifying the smallest modules
  size[k] <- min(front_pareto$len)
  k <- k + 1
  front_pareto %<>%
    filter(len == min(front_pareto$len))
  ## keeping the modules unique
  front_pareto <- front_pareto %>% distinct(za,nest,rnkIndex,len,hits,fp,fn,precision,Recall,F1)
  ## pick up the best module based on the za score.
  index <- which(front_pareto$za == max(front_pareto$za))
  best_module <- front_pareto[index,]
  best_module <- cbind(f,best_module)
  ## Adding the result to the final table.
  final_table <- rbind(final_table,best_module)
}

## Selecting the ami-nsgaII recall of each network
set <- seq(1,nrow(final_table), by = 3)
aminsga2_recall <- NULL; k <- 1
size_final <- NULL
for(j in set){
  aminsga2_recall[k] <- (final_table[j,"hits"] + final_table[j+1,"hits"] + final_table[j+2,"hits"])/30
  size_final[k] <- size[j] + size[j+1] + size[j+2]
  k <- k + 1
}

size_dataset <- length(size_final)
knode_recall <- NULL; bionet_recall <- NULL
## Identifying the Bionet and Knode recall based on Robinson's (2017) results
for(set in 1:size_dataset){
  max <- size_final[set]
  ## Best Knode results 
  knoderesults <- data_sim$overall.knoderesults[,set]
  knoderesults <- order(knoderesults, decreasing = T)[1:max]
  ## Identifying truehists
  truehits <- data_sim$overall.truehits[,set]
  truehits <- which(truehits == 1)
  ## calculating the Knode recall
  tp_knode <- length(intersect(knoderesults,truehits))
  knode_recall[set] <- tp_knode / length(truehits)
  ## Bionet results 
  bionetresults <- data_sim$overall.bionetresults[,set]
  bionetresults <- which(bionetresults == 1)
  ## calculating the Bionet recall
  tp_bionet <- length(intersect(bionetresults,truehits))
  bionet_recall[set] <- tp_bionet / length(truehits)
}

## Identifying the COSINE recall based on Robinson's (2017) results
cosine_results <- NULL
for(set in 1:size_dataset){
  ## update Robinson (2017) data
  max <- size_final[set]
  matrix_adj <- data_sim$overall.adj[[set]]
  matrix_adj <- matrix_adj[[1]]
  matrix_adj <- as.matrix(matrix_adj)
  rownames(matrix_adj) <- seq(1,nrow(matrix_adj))
  colnames(matrix_adj) <- rownames(matrix_adj)
  transcript <- data_sim$overall.pvalues[,set]
  names(transcript) <- rownames(matrix_adj)
  truehits <- data_sim$overall.truehits[,set]
  names(truehits) <- rownames(matrix_adj)
  truehits <- names(which(truehits == 1))
  
  ## Create ppi network ####
  ppi <- graph.adjacency(as.matrix(matrix_adj), weighted = NULL, mode = "directed")
  
  ## Fv = 1 - simulated p-value
  fv <- 1 - as.numeric(transcript)
  
  ## Initialize network weight
  V(ppi)$fc <- rep(0,length(V(ppi)$name))
  
  ## Calculate scaled Nodescore
  scaled_node_score <- NULL
  for(i in 1:length(transcript)){
    fv_i <- fv[i]
    ## Nodescore = Fv - mean(Fv) / stardad deviton(Fv)
    scaled_node_score[i] <- (fv_i - mean(fv)) / sd(fv)
    V(ppi)$fc[i] <- scaled_node_score[i]
  }
  names(scaled_node_score) <- names(transcript)
  
  ## Identifying PPI network
  edges = get.edgelist(ppi)
  e1 <- NULL; e2 <- NULL
  for(i in 1:nrow(edges)){
    x1 <- edges[i,1]
    x2 <- edges[i,2]
    i1 <- which(e2 == x1)
    i2 <- which(e1 == x2)
    test <- intersect(i1,i2)
    if(length(test) == 0){
      e1 <- c(e1,x1)
      e2 <- c(e2,x2)
    }
  }
  PPI <- cbind(as.numeric(e1),as.numeric(e2))
  
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
  klist <- 30; pop_size <- 100
  quantiles <- get_quantiles_PPI(scaled_node_score, scaled_edge_score, PPI, klist, pop_size)
  lambda <- quantiles[[2]]
  
  ## Exectue COSINE for the 5 best lambdas
  best_score <- 0
  for(i in 1:length(lambda)){
    lbd <- lambda[i]
    GA_result <- GA_search_PPI(lambda=lbd,scaled_node_score,scaled_edge_score,PPI,
                               num_iter=300, muCh=0.05, zToR=10, minsize=30)
    
    bsscore <- GA_result$Best_Scores
    if(bsscore > best_score){
      ## Best COSINE results 
      index_min <- order(transcript[GA_result$Subnet], decreasing = F)[1:max]
      subnet <- GA_result$Subnet[index_min]
      tp <- length(intersect(subnet,as.numeric(truehits)))
      fp <- length(setdiff(subnet,as.numeric(truehits)))
      best_score <- bsscore
    }
  }
  
  ## Note the results based in the best lambda value
  fn <- length(truehits) - tp
  recall <- tp/(tp+fn) 
  precision <- tp/(tp+fp)
  f1_score <- 2 * ((precision * recall)/(precision + recall))
  aux <- c(recall,precision,f1_score)
  names(aux) <- c("recall","precision","f1_score")
  cosine_results <- rbind(cosine_results, aux)
}
cossine_recall <- cosine_results$recall

recall_table <- cbind(aminsga2_recall,cossine_recall,knode_recall,bionet_recall)
colnames(recall_table) <- c("AMI-NSGAII","COSSINE","Knode","Bionet")
boxplot(recall_table, main = "50 simulations")

write.csv(recall_table, file = OUTPUT_PATH, row.names = F)
