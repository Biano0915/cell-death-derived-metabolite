library(tidyverse)
# Import csv data
mix_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/Metaboanalyst/Mix/Mix_iden.csv", stringsAsFactors = FALSE)
mix_measure <- mix_measure[, 2:ncol(mix_measure)]

POS_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_POS/3_QI Report/POS_Compound Measurement.csv", stringsAsFactors = FALSE)
POS_measure <- POS_measure[, 14:28]

NEG_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_NEG/3_QI Report/NEG_Compound Measurement.csv", stringsAsFactors = FALSE)
NEG_measure <- NEG_measure[, 14:28]

# checking all the number in the data are numeric
apply(mix_measure, 2, is.numeric)
apply(POS_measure, 2, is.numeric)
apply(NEG_measure, 2, is.numeric)

# checking the number of NA in the data
colSums(is.na(mix_measure))
colSums(is.na(POS_measure))
colSums(is.na(NEG_measure))

# calculate the mean of each group
mix_mean <- mix_measure %>% 
  mutate(DMSO = (mix_measure$DMSO.1 + mix_measure$DMSO.2 + mix_measure$DMSO.3) / 3,
         STS = (mix_measure$STS.1 + mix_measure$STS.2 + mix_measure$STS.3) / 3,
         Fas = (mix_measure$Fas.1 + mix_measure$Fas.2 + mix_measure$Fas.3) / 3,
         RSL3 = (mix_measure$RSL.1 + mix_measure$RSL.2 + mix_measure$RSL.3) / 3,
         Nec = (mix_measure$Nec.1 + mix_measure$Nec.2 + mix_measure$Nec.3) / 3,)
mix_mean <- mix_mean[, 16:20]

pos_mean <- POS_measure %>% 
  mutate(DMSO = (POS_measure$DMSO.1 + POS_measure$DMSO.2 + POS_measure$DMSO.3) / 3,
         STS = (POS_measure$STS.1 + POS_measure$STS.2 + POS_measure$STS.3) / 3,
         Fas = (POS_measure$Fas.1 + POS_measure$Fas.2 + POS_measure$Fas.3) / 3,
         RSL3 = (POS_measure$RSL.1 + POS_measure$RSL.2 + POS_measure$RSL.3) / 3,
         Nec = (POS_measure$Nec.1 + POS_measure$Nec.2 + POS_measure$Nec.3) / 3,)         
pos_mean <- pos_mean[, 16:20]

neg_mean <- NEG_measure %>% 
  mutate(DMSO = (NEG_measure$DMSO.1 + NEG_measure$DMSO.2 + NEG_measure$DMSO.3) / 3,
         STS = (NEG_measure$STS.1 + NEG_measure$STS.2 + NEG_measure$STS.3) / 3,
         Fas = (NEG_measure$Fas.1 + NEG_measure$Fas.2 + NEG_measure$Fas.3) / 3,
         RSL3 = (NEG_measure$RSL.1 + NEG_measure$RSL.2 + NEG_measure$RSL.3) / 3,
         Nec = (NEG_measure$Nec.1 + NEG_measure$Nec.2 + NEG_measure$Nec.3) / 3,)
neg_mean <- neg_mean[, 16:20]

# library
library(FactoMineR)
library(factoextra)
library(ggcorrplot)
library(corrr)



# normalize the data
mix_mean_nor <- scale(mix_mean, center = TRUE, scale = TRUE)


plot(log(mix_mean$RSL3), log(mix_mean$Fas))

# correlation plot
mix_cor <- cor(mix_mean_nor)

ggcorrplot(cor(mix_cor))

