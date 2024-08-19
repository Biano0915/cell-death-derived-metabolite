#----0. library----
library(tidyverse)
library(KEGGREST)

input_path <- "1_Data/20231204 Cell-Death Metabolite_Sup/output/For Pathway/Mix/Unique/"

#----1. import unique ID----
STS_unique_id <- read_csv(paste0(input_path, "STS_mix_unique_ID.csv"))
STS_unique_id <- STS_unique_id[, 2]
colnames(STS_unique_id) <- c("KEGG_ID")

Fas_unique_id <- read_csv(paste0(input_path, "Fas_mix_unique_ID.csv"))
Fas_unique_id <- Fas_unique_id[, 2]
colnames(Fas_unique_id) <- c("KEGG_ID")

RSL_unique_id <- read_csv(paste0(input_path, "RSL_mix_unique_ID.csv"))
RSL_unique_id <- RSL_unique_id[, 2]
colnames(RSL_unique_id) <- c("KEGG_ID")

Nec_unique_id <- read_csv(paste0(input_path, "Nec_mix_unique_ID.csv"))
Nec_unique_id <- Nec_unique_id[, 2]
colnames(Nec_unique_id) <- c("KEGG_ID")

#----2. initialize the list of each group----
STS_list <- list()
Fas_list <- list()
RSL_list <- list()
Nec_list <- list()

test <- function(KEGGID, type = ){
  x <- keggGet(KEGGID)
  if(x[[1]]$NAME != NULL){
    
  }
}