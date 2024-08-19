library(rvest)
library(tidyverse)

URL <- function(CSID){
  url <- paste("http://www.chemspider.com/ibcontent.ashx?csid=", # import the url
               CSID,
               "&type=ds&ds_types=",
               sep = "")
  read <- read_html(url)
  node_1 <- '//*[@id = "DS5"]//a'                                  # give the wanted node
  node_2 <- '//*[@id = "DS1"]//a'
  html_extract_1 <- read %>%                                       # extract the text from the web
    html_nodes(xpath = node_1) %>% 
    html_text()
  html_extract_2 <- read %>% 
    html_nodes(xpath = node_2) %>% 
    html_text()
  html_df_1 <- data.frame(html_extract_1)                        # transfer the data into dataframe format
  html_df_1[] <- lapply(html_df_1[], as.character)
  html_df_2 <- data.frame(html_extract_2)
  html_df_2[] <- lapply(html_df_2[], as.character)
  kegg_row <- which(html_df_1$html_extract_1 == "KEGG")              # find out which row (number) has the string "KEGG"
  hmdb_row <- which(html_df_2$html_extract_2 == "Human Metabolome Database")
  c1 <- grepl("C",                                               # find out whether the KEGG ID is for metabolite, drug, multiple match or doesnt' exist
              substr(html_df_1[kegg_row + 1,], 1, 1)
  )
  c2 <- grepl("C",
              substr(html_df_1[kegg_row + 2,], 1, 1)
  )
  d1 <- grepl("D",
              substr(html_df_1[kegg_row + 1,], 1, 1)
  )
  d2 <- grepl("D",
              substr(html_df_1[kegg_row + 2,], 1, 1)
  )
  h1 <- grepl("HMDB",
              substr(html_df_2[hmdb_row + 1,], 1, 4)
  )
  if(isTRUE(c1 && c2)){
    "MATCH ERROR"
  }
  else if(isTRUE(c1)){
    html_df_1[kegg_row + 1,]
  }
  else if(isTRUE(c2)){
    html_df_1[kegg_row + 2,]
  }
  else if(isTRUE(d1 || d2)){
    "Drug Only"
  }
  else if(isTRUE(h1)){
    html_df_2[hmdb_row + 1,]
  }
  else{
    NA
  }
}

#----1. Converted CSID to KEGG ID----
##----1.1. Negative mode----
b_neg_iden_fifteen_ppm$`Compound ID` <- gsub("CSID", "", b_neg_iden_fifteen_ppm$`Compound ID`)

system.time(
  neg_id <- apply(data.frame(b_neg_iden_fifteen_ppm$`Compound ID`), 1, URL)
)

# merge data
neg_id <- data.frame(b_neg_iden_fifteen_ppm[, 1],
                     KEGG_ID = neg_id,
                     b_neg_iden_fifteen_ppm[, -1])
# export data
write_csv(neg_id, file = "D:/1_Metabolomic/1_Data/20231204 Cell-Death Metabolite_Sup/output/For venn/neg_id_15ppm.csv")

##----1.2. Positive mode----
b_pos_iden_fifteen_ppm$`Compound ID` <- gsub("CSID", "", b_pos_iden_fifteen_ppm$`Compound ID`)


system.time(
  pos_id_1_1000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[1:1000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_1001_2000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[1001:2000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_2001_3000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[2001:3000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_3001_4000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[3001:4000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_4001_5000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[4001:5000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_5001_6000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[5001:6000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_6001_7000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[6001:7000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_7001_8000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[7001:8000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_8001_9000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[8001:9000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_9001_10000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[9001:10000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_10001_11000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[10001:11000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_11001_12000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[11001:12000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_12001_13000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[12001:13000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_13001_14000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[13001:14000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_14001_15000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[14001:15000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_15001_16000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[15001:16000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_16001_17000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[16001:17000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_17001_18000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[17001:18000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_18001_19000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[18001:19000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_19001_20000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[19001:20000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_20001_21000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[20001:21000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_21001_22000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[21001:22000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_22001_23000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[22001:23000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_23001_24000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[23001:24000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_24001_25000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[24001:25000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_25001_26000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[25001:26000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_26001_27000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[26001:27000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_27001_28000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[27001:28000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_28001_29000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[28001:29000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_29001_30000 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[29001:30000]), 1, URL)
)
Sys.sleep(10)

system.time(
  pos_id_30001_30329 <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[30001:30329]), 1, URL)
)

#################################################################################################
#    _________     __     __     __    __    __    __     ______     __    __     ________      #
#   |  _____  \   |  |   |  |   |  \  |  |  |  \  |  |   |_    _|   |  \  |  |   |   __   |     #
#   | |     | |   |  |   |  |   |   \ |  |  |   \ |  |     |  |     |   \ |  |   |  |  |__|     #
#   | |     | |   |  |   |  |   |    \|  |  |    \|  |     |  |     |    \|  |   |  |  ___      #
#   |  -----  |   |  |   |  |   |        |  |        |     |  |     |        |   |  | |_  |     #
#   |   ___  _|   |  |___|  |   |  |\    |  |  |\    |     |  |     |  |\    |   |  |  |  |     #
#   |  |   \ \    |         |   |  | \   |  |  | \   |   __|  |__   |  | \   |   |  |__|  |     #
#   |__|    \_\    \_______/    |__|  \__|  |__|  \__|  |________|  |__|  \__|   |________|     #
#                                                                                               #
#################################################################################################

pos_id <- rbind(data.frame(ID = pos_id_1_1000),
      data.frame(ID = pos_id_1001_2000),
      data.frame(ID = pos_id_2001_3000),
      data.frame(ID = pos_id_3001_4000),
      data.frame(ID = pos_id_4001_5000),
      data.frame(ID = pos_id_5001_6000),
      data.frame(ID = pos_id_6001_7000),
      data.frame(ID = pos_id_7001_8000),
      data.frame(ID = pos_id_8001_9000),
      data.frame(ID = pos_id_9001_10000),
      data.frame(ID = pos_id_10001_11000),
      data.frame(ID = pos_id_11001_12000),
      data.frame(ID = pos_id_12001_13000),
      data.frame(ID = pos_id_13001_14000),
      data.frame(ID = pos_id_14001_15000),
      data.frame(ID = pos_id_15001_16000),
      data.frame(ID = pos_id_16001_17000),
      data.frame(ID = pos_id_17001_18000),
      data.frame(ID = pos_id_18001_19000),
      data.frame(ID = pos_id_19001_20000),
      data.frame(ID = pos_id_20001_21000),
      data.frame(ID = pos_id_21001_22000),
      data.frame(ID = pos_id_22001_23000),
      data.frame(ID = pos_id_23001_24000),
      data.frame(ID = pos_id_24001_25000),
      data.frame(ID = pos_id_25001_26000),
      data.frame(ID = pos_id_26001_27000),
      data.frame(ID = pos_id_27001_28000),
      data.frame(ID = pos_id_28001_29000),
      data.frame(ID = pos_id_29001_30000),
      data.frame(ID = pos_id_30001_30329))
# merge data
pos_id <- data.frame(b_pos_iden_fifteen_ppm[, 1],
                     KEGG_ID = pos_id,
                     b_pos_iden_fifteen_ppm[, -1])

# remove the row contain HMDB in KEGG_ID column
pos_id <- pos_id[!grepl("HMDB", pos_id$ID), ]

# export Data
write_csv(pos_id, file = paste0(output_path, "For venn/pos_id_15ppm.csv"))


