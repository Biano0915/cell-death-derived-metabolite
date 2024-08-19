# This script is for noise reduction of the annotated data of metabolomic
# Threshold, Mass Error <= 5 ppm


#----1. library----
library(tidyverse)
library(rvest)
library(matrixStats) # for median calculation
library(Rdisop)
library(FactoMineR)
library(factoextra)
library(ggcorrplot)

#----2. create the output folder----
dir.create(path = "1_Data/20231204 Cell-Death Metabolite_Sup/output", showWarnings = FALSE)
output_path <- "1_Data/20231204 Cell-Death Metabolite_Sup/output/"

#----3. import raw data----
POS_iden <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_POS/3_QI Report/POS_Compound Identification.csv")
NEG_iden <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_NEG/3_QI Report/NEG_Compound Identification.csv")
POS_measure <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_POS/3_QI Report/POS_Compound Measurement.csv")
NEG_measure <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_NEG/3_QI Report/NEG_Compound Measurement.csv")

#----4. tidy measurement raw data----
col <- c("Compound", 
         "DMSO-1", "DMSO-2", "DMSO-3", 
         "STS-1", "STS-2", "STS-3", 
         "Fas-1", "Fas-2", "Fas-3",
         "RSL-1", "RSL-2", "RSL-3",
         "Nec-1", "Nec-2", "Nec-3",
         "Blank")
POS_measure_clean <- POS_measure[3:9028, c(1, 14:29)]  # extract wanted columns only
colnames(POS_measure_clean) <- col
POS_measure_clean[, 2:17] <- apply(POS_measure_clean[, 2:17], 2, as.numeric)
NEG_measure_clean <- NEG_measure[3:9293, c(1, 14:29)]
colnames(NEG_measure_clean) <- col
NEG_measure_clean[, 2:17] <- apply(NEG_measure_clean[, 2:17], 2, as.numeric)

#----5. calculate mean of each sample & the median of all group----
POS_measure_clean <- POS_measure_clean %>%                              # mean of samples
  mutate(DMSO_mean = rowMeans(POS_measure_clean[, 2:4], na.rm = TRUE),
         STS_mean = rowMeans(POS_measure_clean[, 5:7], na.rm = TRUE),
         Fas_mean = rowMeans(POS_measure_clean[, 8:10], na.rm = TRUE),
         RSL_mean = rowMeans(POS_measure_clean[, 11:13], na.rm = TRUE),
         Nec_mean = rowMeans(POS_measure_clean[, 14:16], na.rm = TRUE))

POS_measure_clean <- POS_measure_clean %>%                              # median of samples
  mutate(All_Sample_Median = rowMedians(as.matrix(POS_measure_clean[, 18:22]), na.rm = TRUE))

NEG_measure_clean <- NEG_measure_clean %>%  
  mutate(DMSO_mean = rowMeans(NEG_measure_clean[, 2:4], na.rm = TRUE),
         STS_mean = rowMeans(NEG_measure_clean[, 5:7], na.rm = TRUE),
         Fas_mean = rowMeans(NEG_measure_clean[, 8:10], na.rm = TRUE),
         RSL_mean = rowMeans(NEG_measure_clean[, 11:13], na.rm = TRUE),
         Nec_mean = rowMeans(NEG_measure_clean[, 14:16], na.rm = TRUE))

NEG_measure_clean <- NEG_measure_clean %>%                              # median of samples
  mutate(All_Sample_Median = rowMedians(as.matrix(NEG_measure_clean[, 18:22]), na.rm = TRUE))

#----6. remove background signal (<= 10x blank)----
POS_measure_rm <- POS_measure_clean[-(which(POS_measure_clean$All_Sample_Median >= 10*POS_measure_clean$Blank)),]
NEG_measure_rm <- NEG_measure_clean[-(which(NEG_measure_clean$All_Sample_Median >= 10*NEG_measure_clean$Blank)),]

#----7. remove background from identification group----
POS_iden_rm <- POS_iden[POS_iden$Compound %in% POS_measure_rm$Compound, ]
NEG_iden_rm <- NEG_iden[NEG_iden$Compound %in% NEG_measure_rm$Compound, ]

#----8. output some data first----
write_csv(POS_iden_rm, 
          paste(output_path, 
                "1_a_POS_identification.csv", 
                sep = ""))
write_csv(NEG_iden_rm, 
          paste(output_path, 
                "1_bNEG_identification.csv", 
                sep = ""))
write_csv(POS_measure_rm, 
          paste(output_path, 
                "2_a_POS_measurement.csv", 
                sep = ""))
write_csv(NEG_measure_rm, 
          paste(output_path, 
                "2_b_NEG_measurement.csv", 
                sep = ""))

#----9. CSID to KEGGID converter----
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

##----9-1. extract CSID only----
POS_ID <- POS_iden_rm$`Compound ID`
POS_CSID <- POS_ID[-(grep("HMDB", POS_ID))] %>% 
  data.frame()
colnames(POS_CSID) <- "CSID"
POS_CSID$CSID <- gsub("CSID", "", as.character(POS_CSID$CSID))   # remove the CSID string

NEG_ID <- NEG_iden_rm$`Compound ID`
NEG_CSID <- NEG_ID[-(grep("HMDB", NEG_ID))] %>% 
  data.frame()
colnames(NEG_CSID) <- "CSID"
NEG_CSID$CSID <- gsub("CSID", "", as.character(NEG_CSID$CSID))   # remove the CSID string

##----9-2. convert CSID----
system.time({                                                    # POS mode ID convert
  x <- apply(data.frame(POS_CSID[1:37192, ]), 1, URL)            # run about 6 hours (20231206 recorded)
})
POSID_converted <- data.frame(x)
POSID_converted <- cbind(POS_CSID, POSID_converted)
colnames(POSID_converted) <- c("Compound ID", "Converted_ID")
POSID_converted$`Compound ID` <- paste("CSID", POSID_converted$`Compound ID`, sep = "")
write_csv(POSID_converted, paste(output_path, 
                                 "3_a_Positive Mode Converted ID.csv", 
                                 sep = ""))

system.time({                                                    # NEG mode ID convert
  y <- apply(data.frame(NEG_CSID[1:6771, ]), 1, URL)             # run about 5037 sec (20231207 recorded)
})
NEGID_converted <- data.frame(y)
NEGID_converted <- cbind(NEG_CSID, NEGID_converted)
colnames(NEGID_converted) <- c("Compound ID", "Converted_ID")
NEGID_converted$`Compound ID` <- paste("CSID", NEGID_converted$`Compound ID`, sep = "")
write_csv(NEGID_converted, paste(output_path, 
                                 "3_b_Negative Mode Converted ID.csv", 
                                 sep = ""))


#----10. remove the row contain HMDBID (col, compound ID)----
POS_iden_CSID <- POS_iden_rm[!(grepl("HMDB", POS_iden_rm$`Compound ID`)), ]
NEG_iden_CSID <- NEG_iden_rm[!(grepl("HMDB", NEG_iden_rm$`Compound ID`)), ]

#----11. merge the converted ID with identification dataframe
POS_iden_CSID <- cbind(POS_iden_CSID, data.frame(POSID_converted$Converted_ID))
colnames(POS_iden_CSID)[18] <- "Converted_ID"
NEG_iden_CSID <- cbind(NEG_iden_CSID, data.frame(NEGID_converted$Converted_ID))
colnames(NEG_iden_CSID)[18] <- "Converted_ID"

#----12. remove the row contain HMDBID & Drug only (col, converted ID)
POS_iden_CSID_clean <- POS_iden_CSID[!(grepl("HMDB", POS_iden_CSID$Converted_ID)), ]
POS_iden_CSID_clean <- POS_iden_CSID_clean[!(grepl("Drug", POS_iden_CSID_clean$Converted_ID)), ]
summary(grepl("M", substr(POS_iden_CSID_clean$Converted_ID, 1, 1)))         # checking how many match error

NEG_iden_CSID_clean <- NEG_iden_CSID[!(grepl("HMDB", NEG_iden_CSID$Converted_ID)), ]
NEG_iden_CSID_clean <- NEG_iden_CSID_clean[!(grepl("Drug", NEG_iden_CSID_clean$Converted_ID)), ]
summary(grepl("M", substr(NEG_iden_CSID_clean$Converted_ID, 1, 1)))         # checking how many match error

#----13. fix match error----
which(POS_iden_CSID_clean$Converted_ID == "MATCH ERROR")        # checking which row has match error
which(NEG_iden_CSID_clean$Converted_ID == "MATCH ERROR")        # checking which row has match error
write_csv(POS_iden_CSID_clean, file = paste(output_path, "4_a_Positive_identification_with_Match_Error.csv", sep = ""))
write_csv(NEG_iden_CSID_clean, file = paste(output_path, "4_b_Negative_identification_with_Match_Error.csv", sep = ""))

#----14. input the fixed file----
tidy_POS_iden <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/output/5_a_Positive_identification_cleaned.csv")
tidy_NEG_iden <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/output/5_b_Negative_identification_cleaned.csv")

#----15. calculate the row has the mass error within +-15ppm----
tidy_NEG_iden$Mass.Error..ppm. <- as.numeric(tidy_NEG_iden$Mass.Error..ppm.)     # make sure all the elements are numeric
tidy_POS_iden$Mass.Error..ppm. <- as.numeric(tidy_POS_iden$Mass.Error..ppm.)

sum(tidy_NEG_iden$Mass.Error..ppm. > -15 & tidy_NEG_iden$Mass.Error..ppm.< 15)   # 4413
sum(tidy_POS_iden$Mass.Error..ppm. > -15 & tidy_POS_iden$Mass.Error..ppm.< 15)   # 4924

sum(tidy_NEG_iden$Mass.Error..ppm. > -10 & tidy_NEG_iden$Mass.Error..ppm.< 10)   # 3318
sum(tidy_POS_iden$Mass.Error..ppm. > -10 & tidy_POS_iden$Mass.Error..ppm.< 10)   # 3678

sum(tidy_NEG_iden$Mass.Error..ppm. > -5 & tidy_NEG_iden$Mass.Error..ppm.< 5)   # 1852
sum(tidy_POS_iden$Mass.Error..ppm. > -5 & tidy_POS_iden$Mass.Error..ppm.< 5)   # 2207

#----16. calculate the molecular weight---- 
weight_count <- function(x){
  getMolecule(x)[[3]]
}
sapply(tidy_POS_iden$Formula, weight_count)






#----17. extract the feature within different mass error range----
five_ppm_POS_iden <- tidy_POS_iden[tidy_POS_iden$Mass.Error..ppm. > -5 & tidy_POS_iden$Mass.Error..ppm.< 5,]
five_ppm_POS_iden <- five_ppm_POS_iden[, !colnames(five_ppm_POS_iden) %in% "Accepted."]   # remove accepted column


ten_ppm_POS_iden <- tidy_POS_iden[tidy_POS_iden$Mass.Error..ppm. > -10 & tidy_POS_iden$Mass.Error..ppm.< 10,]
ten_ppm_POS_iden <- ten_ppm_POS_iden[, !colnames(ten_ppm_POS_iden) %in% "Accepted."]   # remove accepted column


fifteen_ppm_POS_iden <- tidy_POS_iden[tidy_POS_iden$Mass.Error..ppm. > -15 & tidy_POS_iden$Mass.Error..ppm.< 15,]
fifteen_ppm_POS_iden <- fifteen_ppm_POS_iden[, !colnames(fifteen_ppm_POS_iden) %in% "Accepted."]   # remove accepted column


five_ppm_NEG_iden <- tidy_NEG_iden[tidy_NEG_iden$Mass.Error..ppm. > -5 & tidy_NEG_iden$Mass.Error..ppm.< 5,]
five_ppm_NEG_iden <- five_ppm_NEG_iden[, !colnames(five_ppm_NEG_iden) %in% "Accepted."]   # remove accepted column


ten_ppm_NEG_iden <- tidy_NEG_iden[tidy_NEG_iden$Mass.Error..ppm. > -10 & tidy_NEG_iden$Mass.Error..ppm.< 10,]
ten_ppm_NEG_iden <- ten_ppm_NEG_iden[, !colnames(ten_ppm_NEG_iden) %in% "Accepted."]   # remove accepted column


fifteen_ppm_NEG_iden <- tidy_NEG_iden[tidy_NEG_iden$Mass.Error..ppm. > -15 & tidy_NEG_iden$Mass.Error..ppm.< 15,]
fifteen_ppm_NEG_iden <- fifteen_ppm_NEG_iden[, !colnames(fifteen_ppm_NEG_iden) %in% "Accepted."]   # remove accepted column






##----16-2. save the csv files----
write_csv(five_ppm_POS_iden, 
          file = paste(output_path, "6-1_a_5ppm_mass_error_POS_iden.csv", sep = ""))
write_csv(ten_ppm_POS_iden, 
          file = paste(output_path, "6-1_b_10ppm_mass_error_POS_iden.csv", sep = ""))
write_csv(fifteen_ppm_POS_iden, 
          file = paste(output_path, "6-1_c_15ppm_mass_error_POS_iden.csv", sep = ""))
write_csv(five_ppm_NEG_iden, 
          file = paste(output_path, "6-2_a_5ppm_mass_error_NEG_iden.csv", sep = ""))
write_csv(ten_ppm_NEG_iden, 
          file = paste(output_path, "6-2_b_10ppm_mass_error_NEG_iden.csv", sep = ""))
write_csv(fifteen_ppm_NEG_iden, 
          file = paste(output_path, "6-3_a_15ppm_mass_error_NEG_iden.csv", sep = ""))