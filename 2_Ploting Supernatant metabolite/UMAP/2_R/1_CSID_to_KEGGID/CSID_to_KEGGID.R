#----library----
# can install tidyverse to get the rvest package
# still need to library the package seperately
library(rvest)

#----import----
import_data <- read.csv("CSID.csv", header = FALSE)
import_data$V1 <- gsub("CSID", "", as.character(import_data$V1))

#----function----
URL <- function(CSID){
  url <- paste("http://www.chemspider.com/ibcontent.ashx?csid=", # import the url
               CSID,
               "&type=ds&ds_types=",
               sep = "")
  read <- read_html(url)
  node <- '//*[@id = "DS5"]//a'                                  # give the wanted node
  html_extract <- read %>%                                       # extract the text from the web
    html_nodes(xpath = node) %>% 
    html_text()
  html_df <- data.frame(html_extract)                            # transfer the data into dataframe format
  html_df[] <- lapply(html_df[], as.character)
  kegg_row <- which(html_df$html_extract == "KEGG")              # find out which row (number) has the string "KEGG"
  c1 <- grepl("C",                                               # find out whether the KEGG ID is for metabolite, drug, multiple match or doesnt' exist
              substr(html_df[kegg_row + 1,], 1, 1)
              )
  c2 <- grepl("C",
              substr(html_df[kegg_row + 2,], 1, 1)
              )
  d1 <- grepl("D",
              substr(html_df[kegg_row + 1,], 1, 1)
              )
  d2 <- grepl("D",
              substr(html_df[kegg_row + 2,], 1, 1)
              )
  if(c1 && c2 == TRUE){
    "MATCH ERROR"
  }
  else if(c1 == TRUE){
    html_df[kegg_row + 1,]
  }
  else if(c2 == TRUE){
    html_df[kegg_row + 2,]
  }
  else if(d1 || d2){
    "Drug Only"
  }
  else{
    NA
  }
}
#----Convert----
kegg_id_1_100 <- apply(data.frame(import_data[1:100,]), 1, URL)
kegg_id_101_200 <- apply(data.frame(import_data[101:200,]), 1, URL)
kegg_id_201_300 <- apply(data.frame(import_data[201:300,]), 1, URL)
kegg_id_301_400 <- apply(data.frame(import_data[301:400,]), 1, URL)
kegg_id_401_500 <- apply(data.frame(import_data[401:500,]), 1, URL)
kegg_id_501_600 <- apply(data.frame(import_data[501:600,]), 1, URL)
kegg_id_601_700 <- apply(data.frame(import_data[601:700,]), 1, URL)
kegg_id_701_756 <- apply(data.frame(import_data[701:756,]), 1, URL)

kegg_id <- rbind(
  data.frame(kegg_id = kegg_id_1_100, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_101_200, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_201_300, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_301_400, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_401_500, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_501_600, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_601_700, stringsAsFactors = FALSE),
  data.frame(kegg_id = kegg_id_701_756, stringsAsFactors = FALSE)
  )
write.csv(kegg_id, file = "transfered kegg id")


