# This script is use for scrape the necessary from HMDB website #

#----1. Load the necessary packages----
library(rvest)
library(tidyverse)

#----2. read the data----
neg_iden <- read_csv("20240131 Ploting/Processed/3_neg_HMDB_iden_merged.csv")
pos_iden <- read_csv("20240131 Ploting/Processed/3_pos_HMDB_iden_merged.csv")

#----3. function for scrape the data----
URL <- function(HMDB_Link){
  # read the html
  read <- read_html(HMDB_Link)
  
  # setup nodes
  Common_name_node <- '/html/body/main/table/tbody[1]/tr[9]/td/strong'    # common name
  Class_node <- '/html/body/main/table/tbody[1]/tr[26]/td'                # class
  Subclass_node <- '/html/body/main/table/tbody[1]/tr[27]/td'             # subclass
  Status_node <- '/html/body/main/table/tbody[1]/tr[3]/td'                # status
  Direct_parent_node <- '/html/body/main/table/tbody[1]/tr[28]/td'        # direct parent
  
  # extract the text
  common_name <- read %>%
    html_nodes(xpath = Common_name_node) %>%
    html_text() %>%
    data.frame("common_name" = .)
  class <- read %>%
    html_nodes(xpath = Class_node) %>%
    html_text() %>%
    data.frame("class" = .)
  subclass <- read %>%
    html_nodes(xpath = Subclass_node) %>%
    html_text() %>%
    data.frame("subclass" = .)
  status <- read %>%
    html_nodes(xpath = Status_node) %>%
    html_text() %>%
    data.frame("status" = .)
  direct_parent <- read %>%
    html_nodes(xpath = Direct_parent_node) %>%
    html_text() %>%
    data.frame("direct_parent" = .)
  HMDBID <- data.frame("HMDB_ID" = str_extract(HMDB_Link, "HMDB\\d+"))
  
  # combine the HMDBID, name, class, subclass
  profile <- cbind(HMDBID, status, common_name, direct_parent, class, subclass)
  return(profile) # return the result
}

#----4. scrape the data----
##----4-1. negative data----
neg_result <- list()                                         # Initialize the list
neg_error <- c()                                             # Initialize the error list

for(i in 1 : nrow(neg_iden)){
  # print the progress
  cat("\n----", "\nProcessing (", i, "/", nrow(neg_iden), ")...\n")
  # scrape the data
  tryCatch(                                     # prevent error turn off the loop
    {
      neg_result[[i]] <- URL(neg_iden$Link[i])    # scrape the data then store in the list
    }, 
    error = function(e)                         # if error occur
    {
      cat("\nERROR OCCUR: Num.", i, ", input the data with NA\n")
      neg_error <<- c(neg_error, as.numeric(i))              # store the error number
      neg_result[[i]] <<- data.frame("HMDB_ID" = str_extract(neg_iden$Link[i], "HMDB\\d+"),
                                     "status" = "error",
                                     "common_name" = NA,
                                     "direct_parent" = NA,
                                     "class" = NA, 
                                     "subclass" = NA,
                                     "physiological_effect" = NA)
    }
  )
}

# error re-run
for(i in neg_error){
  # print the progress
  cat("\n----", "\nProcessing (", i, "/", nrow(neg_iden), ")...\n")
  # scrape the data
  tryCatch(                                     # prevent error turn off the loop
    {
      neg_result[[i]] <- URL(neg_iden$Link[i])    # scrape the data then store in the list
      # remove from the error ID
      neg_error <- neg_error[-which(neg_error == i)]
    }, 
    error = function(e)                         # if error occur
    {
      cat("\nERROR OCCUR: Num.", i, ", input the data with NA\n")
      neg_result[[i]] <<- data.frame("HMDB_ID" = str_extract(neg_iden$Link[i], "HMDB\\d+"),
                                     "status" = "error",
                                     "common_name" = NA,
                                     "direct_parent" = NA,
                                     "class" = NA, 
                                     "subclass" = NA,
                                     "physiological_effect" = NA)
    }
  )
}




neg_result <- do.call(rbind, neg_result)          # combine the result

write_csv(neg_result, file = "20240131 Ploting/Processed/3-1_neg_HMDB_parent.csv") # save the data as csv file

##----4-2. positive data----
pos_result <- list()                                         # Initialize the list
pos_error <- c()                                             # Initialize the error list

for(i in 1 : nrow(pos_iden)){
  # print the progress
  cat("\n----", "\nProcessing (", i, "/", nrow(pos_iden), ")...\n")
  # scrape the data
  tryCatch(                                     # prevent error turn off the loop
    {
      pos_result[[i]] <- URL(pos_iden$Link[i])    # scrape the data then store in the list
    }, 
    error = function(e)                         # if error occur
    {
      cat("\nERROR OCCUR: Num.", i, ", input the data with NA\n")
      pos_error <<- c(pos_error, as.numeric(i))              # store the error number
      pos_result[[i]] <<- data.frame("HMDB_ID" = str_extract(pos_iden$Link[i], "HMDB\\d+"), 
                                     "status" = "error",
                                     "common_name" = NA,
                                     "direct_parent" = NA,
                                     "class" = NA, 
                                     "subclass" = NA,
                                     "physiological_effect" = NA)
    }
  )
}

# error re-run
for(i in pos_error){
  # print the progress
  cat("\n----", "\nProcessing (", i, "/", nrow(pos_iden), ")...\n")
  # scrape the data
  tryCatch(                                     # prevent error turn off the loop
    {
      pos_result[[i]] <- URL(pos_iden$Link[i])    # scrape the data then store in the list
      # remove from the error ID
      pos_error <- pos_error[-which(pos_error == i)]
    }, 
    error = function(e)                         # if error occur
    {
      cat("\nERROR OCCUR: Num.", i, ", input the data with NA\n")
      pos_result[[i]] <<- data.frame("HMDB_ID" = str_extract(pos_iden$Link[i], "HMDB\\d+"), 
                                     "status" = "error",
                                     "common_name" = NA,
                                     "direct_parent" = NA,
                                     "class" = NA, 
                                     "subclass" = NA,
                                     "physiological_effect" = NA)
    }
  )
}

pos_result <- do.call(rbind, pos_result)          # combine the result

write_csv(pos_result, file = "20240131 Ploting/Processed/3-1_pos_HMDB_parent.csv") # save the data as csv file
