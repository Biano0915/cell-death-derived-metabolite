# This script is for noise reduction of the annotated data of metabolomic
# Threshold 1, Mass Error <= 5 ppm & >= -5 ppm
# Threshold 2, Mass Error <= 15 ppm & >= -15 ppm


#----1. library----
library(tidyverse)
library(rvest)
library(matrixStats) # for median calculation
library(Rdisop)

#----2. create the output folder----
dir.create(path = "1_Data/20231204 Cell-Death Metabolite_Sup/output", showWarnings = FALSE)
output_path <- "1_Data/20231204 Cell-Death Metabolite_Sup/output/"

#----3. import the identified data----
a_pos_iden_rm <- read_csv(paste(output_path, "1_a_POS_identification.csv", sep = ""))

a_neg_iden_rm <- read_csv(paste(output_path, "1_b_NEG_identification.csv", sep = ""))

#----4. Remove the row contain HMDB ID then filter the row with Mass error <=15 & >=-15----
##----4.1. Positive Mode----
b_pos_iden_fifteen_ppm <- a_pos_iden_rm %>% 
  filter(!grepl("HMDB", `Compound ID`)) %>% 
  filter(`Mass Error (ppm)` <= 15 & `Mass Error (ppm)` >= -15)

##---4.2. negative Mode----
b_neg_iden_fifteen_ppm <- a_neg_iden_rm %>%
  filter(!grepl("HMDB", `Compound ID`)) %>% 
  filter(`Mass Error (ppm)` <= 15 & `Mass Error (ppm)` >= -15)

#----5. Convert CSID to KEGG ID (only run once, need around 6~7 hours)----
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
##----5.1. Negative mode----
b_neg_iden_fifteen_ppm$`Compound ID` <- gsub("CSID", "", b_neg_iden_fifteen_ppm$`Compound ID`)

system.time(
  neg_id <- apply(data.frame(b_neg_iden_fifteen_ppm$`Compound ID`), 1, URL)
)

# merge data
neg_id <- data.frame(b_neg_iden_fifteen_ppm[, 1],
                     KEGG_ID = neg_id,
                     b_neg_iden_fifteen_ppm[, -1])

# export data
write_csv(neg_id, file = paste0(output_path, "For venn/neg_id_15ppm.csv"))

##----5.2. Positive mode----
b_pos_iden_fifteen_ppm$`Compound ID` <- gsub("CSID", "", b_pos_iden_fifteen_ppm$`Compound ID`)

# setup the loop
total_compounds <- nrow(b_pos_iden_fifteen_ppm)
batch_size <- 1000
pos_id <- data.frame(ID = character())

# loop
for (i in seq(1, total_compounds, by = batch_size)) {
  start <- i
  end <- min(i + batch_size - 1, total_compounds)
  
  # print the progress
  cat("Processing", start, "to", end, "of", total_compounds, "...\n")
  
  # measure the time
  start_time <- Sys.time()
  
  # Apply the URL function to the current batch and append to pos_id
  batch_result <- apply(data.frame(b_pos_iden_fifteen_ppm$`Compound ID`[start:end]), 1, URL)
  
  pos_id <- rbind(pos_id, data.frame(ID = batch_result))
  
  end_time <- Sys.time()
  cat("Time taken: ", end_time - start_time, "seconds\n")
  
  # Sleep for 10 seconds
  Sys.sleep(10)
}

##----5.3 import data (this step is for the second time running the code)----
pos_id <- read_csv(paste0(output_path, "For venn/1_Converted/pos_id_15ppm.csv"))
neg_id <- read_csv(paste0(output_path, "For venn/1_Converted/neg_id_15ppm.csv"))
#----6. Import the measurement data-----
c_pos_measure <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_POS/3_QI Report/POS_Compound Measurement.csv")
c_neg_measure <- read_csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_NEG/3_QI Report/NEG_Compound Measurement.csv")

#----7. tidy the measurement file----
col <- c("Compound", 
         "DMSO-1", "DMSO-2", "DMSO-3", 
         "STS-1", "STS-2", "STS-3", 
         "Fas-1", "Fas-2", "Fas-3",
         "RSL-1", "RSL-2", "RSL-3",
         "Nec-1", "Nec-2", "Nec-3",
         "Blank")
##----7.1 Positive mode----
c_pos_measure_clean <- c_pos_measure[3:9028, c(1, 14:29)]  # extract wanted columns only
colnames(c_pos_measure_clean) <- col  # rename the columns
c_pos_measure_clean[, 2:17] <- apply(c_pos_measure_clean[, 2:17], 2, as.numeric) # make sure all the element are numeric

##----7.2 Negative mode----
c_neg_measure_clean <- c_neg_measure[3:9293, c(1, 14:29)]  # extract wanted columns only
colnames(c_neg_measure_clean) <- col  # rename the columns
c_neg_measure_clean[, 2:17] <- apply(c_neg_measure_clean[, 2:17], 2, as.numeric) # make sure all the element are numeric

#----8. Remove the background signal----
##----8.1 Calculate the median of the all samples----
###----8.1.1 Positive Mode----
c_pos_measure_clean <- c_pos_measure_clean %>%                              # mean of samples
  mutate(DMSO_mean = rowMeans(c_pos_measure_clean[, 2:4], na.rm = TRUE),
         STS_mean = rowMeans(c_pos_measure_clean[, 5:7], na.rm = TRUE),
         Fas_mean = rowMeans(c_pos_measure_clean[, 8:10], na.rm = TRUE),
         RSL_mean = rowMeans(c_pos_measure_clean[, 11:13], na.rm = TRUE),
         Nec_mean = rowMeans(c_pos_measure_clean[, 14:16], na.rm = TRUE))

c_pos_measure_clean <- c_pos_measure_clean %>%                              # median of samples
  mutate(All_Sample_Median = rowMedians(as.matrix(c_pos_measure_clean[, 18:22]), na.rm = TRUE))

###----8.1.2 Negative Mode----
c_neg_measure_clean <- c_neg_measure_clean %>%                              # mean of samples
  mutate(DMSO_mean = rowMeans(c_neg_measure_clean[, 2:4], na.rm = TRUE),
         STS_mean = rowMeans(c_neg_measure_clean[, 5:7], na.rm = TRUE),
         Fas_mean = rowMeans(c_neg_measure_clean[, 8:10], na.rm = TRUE),
         RSL_mean = rowMeans(c_neg_measure_clean[, 11:13], na.rm = TRUE),
         Nec_mean = rowMeans(c_neg_measure_clean[, 14:16], na.rm = TRUE))

c_neg_measure_clean <- c_neg_measure_clean %>%                              # median of samples
  mutate(All_Sample_Median = rowMedians(as.matrix(c_neg_measure_clean[, 18:22]), na.rm = TRUE))

##----8.2 Remove background signal----
###----8.2.1 Positive Mode----
c_pos_measure_rm <- c_pos_measure_clean[-(which(c_pos_measure_clean$All_Sample_Median >= 10*c_pos_measure_clean$Blank)),]

###----8.2.2 Negative Mode----
c_neg_measure_rm <- c_neg_measure_clean[-(which(c_neg_measure_clean$All_Sample_Median >= 10*c_neg_measure_clean$Blank)),]

#----9. UMAP----
library(umap)
plot.umap <- function(x, labels,
                      main= "UMAP",
                      colors=c("#ff7f00", "#ff7f00", "#ff7f00", 
                               "#e377c2", "#e377c2", "#e377c2", 
                               "#17becf", "#17becf", "#17becf", 
                               "#bcbd22", "#bcbd22", "#bcbd22", 
                               "#1f77b4", "#1f77b4", "#1f77b4"),
                      pad=0.5, cex=1, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1, cex.legend=0.85) {
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  } 
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

##----9.1 Positive mode----
set.seed(1)
pos_umap_data <- c_pos_measure_rm[, 2:16] %>% 
  na.omit() %>% 
  t() %>% 
  data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>% 
  umap()
pos_umap_labels <- as.factor(colnames(c_pos_measure_rm[2:16]))

# plot positive mode UMAP
plot.umap(pos_umap_data, pos_umap_labels, main = "UMAP of Positive Mode")

##----9.2 Negative mode----
set.seed(1)
neg_umap_data <- c_neg_measure_rm[, 2:16] %>% 
  na.omit() %>% 
  t() %>% 
  data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>% 
  umap()
neg_umap_labels <- as.factor(colnames(c_neg_measure_rm[2:16]))

# plot negative mode UMAP
plot.umap(neg_umap_data, neg_umap_labels, main = "UMAP of Negative Mode")

##----9.3 Mix mode----
###----9.3.1 combine postive and negative mode data----
mix_measure <- rbind(c_pos_measure_rm, c_neg_measure_rm)
set.seed(1)
mix_measure_umap_data <- mix_measure[, 2:16] %>% 
  na.omit() %>% 
  t() %>% 
  data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>% 
  umap()
mix_measure_umap_labels <- as.factor(colnames(mix_measure[2:16]))

# plot mix data UMAP
plot.umap(mix_measure_umap_data, mix_measure_umap_labels, main = "UMAP of Mix Data")

#----10. Filter the identify data (remove drug and match error)----
##----10.1. Positive mode----
# remove Drug row
pos_id <- pos_id[-grep("Drug", pos_id$ID), ]

# make sure all the element in mass error column are numeric
pos_id$Mass.Error..ppm. <- as.numeric(pos_id$Mass.Error..ppm.)

# remove the duplicate rows (has the same compound and ID), keep the one with smallest mass error
# if the mass error are the same, keep one row randomly
pos_id <- pos_id %>% 
  group_by(Compound, ID) %>% 
  filter(abs(Mass.Error..ppm.) == min(abs(Mass.Error..ppm.))) %>%
  slice_sample(n = 1)

# save the processed data, for Match Error check (manually)
write.csv(pos_id, file = paste0(output_path, "For venn/2_Without Drug/pos_id_15ppm.csv"), row.names = FALSE)

##----10.2. Negative mode----
# remove Drug row
neg_id <- neg_id[-grep("Drug", neg_id$KEGG_ID), ]

# make sure all the element in mass error column are numeric
neg_id$Mass.Error..ppm. <- as.numeric(neg_id$Mass.Error..ppm.)

# remove the duplicate rows (has the same compound and ID), keep the one with smallest mass error
# if the mass error are the same, keep one row randomly
neg_id <- neg_id %>% 
  group_by(Compound, KEGG_ID) %>% 
  filter(abs(Mass.Error..ppm.) == min(abs(Mass.Error..ppm.))) %>%
  slice_sample(n = 1)

# save the processed data, for Match Error check (manually)
write.csv(neg_id, file = paste0(output_path, "For venn/2_Without Drug/neg_id_15ppm.csv"), row.names = FALSE)

######################################################
#                                                    #
#                                                    #
#     Stop here, correct the match error manually    #
#                                                    #
#                                                    #
######################################################

#----11. Import the corrected match error data----
pos_id <- read_csv(paste0(output_path, "For venn/3_Without Match Error/pos_id_15ppm.csv"))
neg_id <- read.csv(paste0(output_path, "For venn/3_Without Match Error/neg_id_15ppm.csv"))

#----12. Remove the same annotation
##----12.1. Positive mode----
pos_id <- pos_id %>% 
  group_by(Compound, ID) %>% 
  filter(abs(Mass.Error..ppm.) == min(abs(Mass.Error..ppm.))) %>%
  slice_sample(n = 1)

##----12.2. Negative mode----
neg_id <- neg_id %>% 
  group_by(Compound, KEGG_ID) %>% 
  filter(abs(Mass.Error..ppm.) == min(abs(Mass.Error..ppm.))) %>%
  slice_sample(n = 1)

#----13. Summarize the annotation data----
##----13.1. Positive mode----
pos_id_sum <- pos_id %>% 
  group_by(Compound) %>%
  summarise(ID = paste(ID, collapse = ", "))

##----13.2. Negative mode----
neg_id_sum <- neg_id %>% 
  group_by(Compound) %>%
  summarise(ID = paste(KEGG_ID, collapse = ", "))

#----14. Remove the measure data with no annotation----
##----14.1. Positive mode----
pos_measure <- merge(c_pos_measure_rm, pos_id_sum, by = "Compound")

##----14.2. Negative mode----
neg_measure <- merge(c_neg_measure_rm, neg_id_sum, by = "Compound")

#----15. Calculate the feature with more than 2x fold change compare to DMSO----
##----15.1. Positive mode----
###----15.1.1. STS----
STS_pos_measure <- pos_measure[, c("Compound", "DMSO_mean", "STS_mean", "ID")]
STS_pos_measure$fold_change <- STS_pos_measure$STS_mean / STS_pos_measure$DMSO_mean
STS_pos_measure <- STS_pos_measure[c(STS_pos_measure$fold_change >= 2), ]
STS_pos_ID <- STS_pos_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.1.2. Fas----
Fas_pos_measure <- pos_measure[, c("Compound", "DMSO_mean", "Fas_mean", "ID")]
Fas_pos_measure$fold_change <- Fas_pos_measure$Fas_mean / Fas_pos_measure$DMSO_mean
Fas_pos_measure <- Fas_pos_measure[c(Fas_pos_measure$fold_change >= 2), ]
Fas_pos_ID <- Fas_pos_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.1.3. RSL----
RSL_pos_measure <- pos_measure[, c("Compound", "DMSO_mean", "RSL_mean", "ID")]
RSL_pos_measure$fold_change <- RSL_pos_measure$RSL_mean / RSL_pos_measure$DMSO_mean
RSL_pos_measure <- RSL_pos_measure[c(RSL_pos_measure$fold_change >= 2), ]
RSL_pos_ID <- RSL_pos_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.1.4. Nec----
Nec_pos_measure <- pos_measure[, c("Compound", "DMSO_mean", "Nec_mean", "ID")]
Nec_pos_measure$fold_change <- Nec_pos_measure$Nec_mean / Nec_pos_measure$DMSO_mean
Nec_pos_measure <- Nec_pos_measure[c(Nec_pos_measure$fold_change >= 2), ]
Nec_pos_ID <- Nec_pos_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()


###----15.1.5. Save the data----
write.csv(STS_pos_ID, file = paste0(output_path, "For Pathway/POS/STS_pos_ID.csv"), row.names = FALSE)
write.csv(Fas_pos_ID, file = paste0(output_path, "For Pathway/POS/Fas_pos_ID.csv"), row.names = FALSE)
write.csv(RSL_pos_ID, file = paste0(output_path, "For Pathway/POS/RSL_pos_ID.csv"), row.names = FALSE)
write.csv(Nec_pos_ID, file = paste0(output_path, "For Pathway/POS/Nec_pos_ID.csv"), row.names = FALSE)

##----15.2. Negative mode----
###----15.2.1. STS----
STS_neg_measure <- neg_measure[, c("Compound", "DMSO_mean", "STS_mean", "ID")]
STS_neg_measure$fold_change <- STS_neg_measure$STS_mean / STS_neg_measure$DMSO_mean
STS_neg_measure <- STS_neg_measure[c(STS_neg_measure$fold_change >= 2), ]
STS_neg_ID <- STS_neg_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.2.2. Fas----
Fas_neg_measure <- neg_measure[, c("Compound", "DMSO_mean", "Fas_mean", "ID")]
Fas_neg_measure$fold_change <- Fas_neg_measure$Fas_mean / Fas_neg_measure$DMSO_mean
Fas_neg_measure <- Fas_neg_measure[c(Fas_neg_measure$fold_change >= 2), ]
Fas_neg_ID <- Fas_neg_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.2.3. RSL----
RSL_neg_measure <- neg_measure[, c("Compound", "DMSO_mean", "RSL_mean", "ID")]
RSL_neg_measure$fold_change <- RSL_neg_measure$RSL_mean / RSL_neg_measure$DMSO_mean
RSL_neg_measure <- RSL_neg_measure[c(RSL_neg_measure$fold_change >= 2), ]
RSL_neg_ID <- RSL_neg_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.2.4. Nec----
Nec_neg_measure <- neg_measure[, c("Compound", "DMSO_mean", "Nec_mean", "ID")]
Nec_neg_measure$fold_change <- Nec_neg_measure$Nec_mean / Nec_neg_measure$DMSO_mean
Nec_neg_measure <- Nec_neg_measure[c(Nec_neg_measure$fold_change >= 2), ]
Nec_neg_ID <- Nec_neg_measure$ID %>% 
  strsplit(split = ", ") %>% 
  unlist() %>% 
  unique()

###----15.2.5. Save the data----
write.csv(STS_neg_ID, file = paste0(output_path, "For Pathway/NEG/STS_neg_ID.csv"), row.names = FALSE)
write.csv(Fas_neg_ID, file = paste0(output_path, "For Pathway/NEG/Fas_neg_ID.csv"), row.names = FALSE)
write.csv(RSL_neg_ID, file = paste0(output_path, "For Pathway/NEG/RSL_neg_ID.csv"), row.names = FALSE)
write.csv(Nec_neg_ID, file = paste0(output_path, "For Pathway/NEG/Nec_neg_ID.csv"), row.names = FALSE)

#----16. Feature Venn Diagram----
library(VennDiagram)

##----16.1 Positive mode----
###----16.1.1. set.seed----
set.seed(12)
color = c("#1f77b4","#e377c2", "#bcbd22", "#17becf")
InApop_P <- STS_pos_measure$Compound %>% 
  na.omit()
ExApop_P <- Fas_pos_measure$Compound %>% 
  na.omit()
Fer_P <- RSL_pos_measure$Compound %>% 
  na.omit()
Nec_P <- Nec_pos_measure$Compound %>% 
  na.omit()

###----16.1.2. Venn Diagram----
ven_pos <- list(InApop_P, ExApop_P, Fer_P, Nec_P)

venn_plot_pos <- venn.diagram(x = ven_pos,
             category.names = c("InApo", "ExApo", "Fer", "Nec"),
             filename = paste0(output_path, "For venn/4_Result//Feature/Positive_Venn.png"),
             output = TRUE,
             imagetype = "png",
             scaled = TRUE,
             col = "black",
             fill = color,
             cat.col = color,
             cat.cex = 2,
             margin = 0.15,
             cex = 1.5,
             cat.fontfamily = "sans")

library(png)
img <- readPNG("Positive_Venn.png")
grid.newpage()
grid.raster(img)

##----16.2 Negative mode----
###----16.2.1. set.seed----
set.seed(12)
InApop_N <- STS_neg_measure$Compound %>%
  na.omit()
ExApop_N <- Fas_neg_measure$Compound %>%
  na.omit()
Fer_N <- RSL_neg_measure$Compound %>%
  na.omit()
Nec_N <- Nec_neg_measure$Compound %>%
  na.omit()

###----16.2.2. Venn Diagram----
ven_neg <- list(InApop_N, ExApop_N, Fer_N, Nec_N)

venn_plot_neg <- venn.diagram(x = ven_neg,
             category.names = c("InApo", "ExApo", "Fer", "Nec"),
             filename = paste0(output_path, "For venn/4_Result//Feature/Negative_Venn.png"),
             output = TRUE,
             imagetype = "png",
             scaled = TRUE,
             col = "black",
             fill = color,
             cat.col = color,
             cat.cex = 2,
             margin = 0.15,
             cat.fontfamily = "sans",
             cex = 1.5)

img <- readPNG("Negative_Venn.png")
grid.newpage()
grid.raster(img)


#----17. Annotation Venn Diagram----
library(VennDiagram)

##----17.1 Positive mode----
set.seed(12)
Nec_pos_ID <- Nec_pos_ID %>% 
  na.omit()
ven_pos_annotated <- list(STS_pos_ID, Fas_pos_ID, RSL_pos_ID, Nec_pos_ID)

venn_plot_pos_annotated <- venn.diagram(x = ven_pos_annotated,
             category.names = c("InApo", "ExApo", "Fer", "Nec"),
             filename = paste0(output_path, "For venn/4_Result/Annotated/Positive_Venn.png"),
             output = TRUE,
             imagetype = "png",
             scaled = TRUE,
             col = "black",
             fill = color,
             cat.col = color,
             cat.cex = 2,
             margin = 0.15,
             cat.fontfamily = "sans")

##----17.2 Negative mode----
set.seed(12)
STS_neg_ID <- STS_neg_ID %>%
  na.omit()
Fas_neg_ID <- Fas_neg_ID %>%
  na.omit()
RSL_neg_ID <- RSL_neg_ID %>%
  na.omit()
Nec_neg_ID <- Nec_neg_ID %>%
  na.omit()
ven_neg_annotated <- list(STS_neg_ID, Fas_neg_ID, RSL_neg_ID, Nec_neg_ID)

venn_plot_neg_annotated <- venn.diagram(x = ven_neg_annotated,
             category.names = c("InApo", "ExApo", "Fer", "Nec"),
             filename = paste0(output_path, "For venn/4_Result/Annotated/Negative_Venn.png"),
             output = TRUE,
             imagetype = "png",
             scaled = TRUE,
             col = "black",
             fill = color,
             cat.col = color,
             cat.cex = 2,
             margin = 0.15,
             cat.fontfamily = "sans")

##----17.3 Mix mode----
set.seed(12)
STS_mix_ID <- c(STS_pos_ID, STS_neg_ID) %>%
  na.omit() %>% 
  unique()
Fas_mix_ID <- c(Fas_pos_ID, Fas_neg_ID) %>%
  na.omit() %>% 
  unique()
RSL_mix_ID <- c(RSL_pos_ID, RSL_neg_ID) %>%
  na.omit() %>% 
  unique()
Nec_mix_ID <- c(Nec_pos_ID, Nec_neg_ID) %>%
  na.omit() %>% 
  unique()
ven_mix_annotated <- list(STS_mix_ID, Fas_mix_ID, RSL_mix_ID, Nec_mix_ID)

venn_plot_mix_annotated <- venn.diagram(x = ven_mix_annotated,
             category.names = c("InApo", "ExApo", "Fer", "Nec"),
             filename = paste0(output_path, "For venn/4_Result/Annotated/Mix_Venn.png"),
             output = TRUE,
             imagetype = "png",
             scaled = TRUE,
             col = "black",
             fill = color,
             cat.col = color,
             cat.cex = 2,
             margin = 0.15,
             cat.fontfamily = "sans",
             cex = 1.5)

###----17.3.1 save mix mode data----
write.csv(STS_mix_ID, paste0(output_path, "For Pathway/Mix/Total/STS_mix_ID.csv"))
write.csv(Fas_mix_ID, paste0(output_path, "For Pathway/Mix/Total/Fas_mix_ID.csv"))
write.csv(RSL_mix_ID, paste0(output_path, "For Pathway/Mix/Total/RSL_mix_ID.csv"))
write.csv(Nec_mix_ID, paste0(output_path, "For Pathway/Mix/Total/Nec_mix_ID.csv"))

#----18. Get the unique ID for each group----
STS_mix_unique_ID <- STS_mix_ID %>% 
  setdiff(c(Fas_mix_ID, RSL_mix_ID, Nec_mix_ID))
Fas_mix_unique_ID <- Fas_mix_ID %>%
  setdiff(c(STS_mix_ID, RSL_mix_ID, Nec_mix_ID))
RSL_mix_unique_ID <- RSL_mix_ID %>%
  setdiff(c(STS_mix_ID, Fas_mix_ID, Nec_mix_ID))
Nec_mix_unique_ID <- Nec_mix_ID %>%
  setdiff(c(STS_mix_ID, Fas_mix_ID, RSL_mix_ID))

##----18.1 save unique ID----
write.csv(STS_mix_unique_ID, paste0(output_path, "For Pathway/Mix/Unique/STS_mix_unique_ID.csv")) #54
write.csv(Fas_mix_unique_ID, paste0(output_path, "For Pathway/Mix/Unique/Fas_mix_unique_ID.csv")) #55
write.csv(RSL_mix_unique_ID, paste0(output_path, "For Pathway/Mix/Unique/RSL_mix_unique_ID.csv")) #249
write.csv(Nec_mix_unique_ID, paste0(output_path, "For Pathway/Mix/Unique/Nec_mix_unique_ID.csv")) #109
