library(tidyverse)
library(corrplot)
library(ggfortify)
library(factoextra)
library(ggforce)

# Load the data
neg_measurement <- read_csv("0_raw_data/NEG HMDB output/NEG HMDB Measurement.csv")
pos_measurement <- read_csv("0_raw_data/POS HMDB output/POS HMDB Measurement.csv")

# remove the unwanted columns (2~6 & 22)
neg_measurement <- neg_measurement %>% select(-c(2:6, 22))
pos_measurement <- pos_measurement %>% select(-c(2:6, 22))

# rename the columns then remove the first two rows
colnames(neg_measurement) <- c("feature", 
                               "Con.1", "Con.2", "Con.3", 
                               "InApop.1", "InApop.2", "InApop.3",
                               "ExApop.1", "ExApop.2", "ExApop.3",
                               "Fer.1", "Fer.2", "Fer.3",
                               "Nec.1", "Nec.2", "Nec.3")
colnames(pos_measurement) <- c("feature", 
                               "Con.1", "Con.2", "Con.3", 
                               "InApop.1", "InApop.2", "InApop.3",
                               "ExApop.1", "ExApop.2", "ExApop.3",
                               "Fer.1", "Fer.2", "Fer.3",
                               "Nec.1", "Nec.2", "Nec.3")
neg_measurement <- neg_measurement[-c(1, 2), ]
pos_measurement <- pos_measurement[-c(1, 2), ]

# transfer the element into numeric
neg_measurement[, 2:16] <- sapply(neg_measurement[, 2:16], as.numeric)
pos_measurement[, 2:16] <- sapply(pos_measurement[, 2:16], as.numeric)

# change tibble to dataframe + change row names with feature, then remove the feature column
neg_measurement <- as.data.frame(neg_measurement)
rownames(neg_measurement) <- neg_measurement$feature
neg_measurement <- neg_measurement[, -1]

pos_measurement <- as.data.frame(pos_measurement)
rownames(pos_measurement) <- pos_measurement$feature
pos_measurement <- pos_measurement[, -1]

# transpose the dataframe
neg_measurement <- t(neg_measurement) %>% 
  as.data.frame()
pos_measurement <- t(pos_measurement) %>% 
  as.data.frame()

# generate the mix dataframe
mix_measurement <- cbind(neg_measurement, pos_measurement)

# add group column
neg_measurement$Group <- gsub("\\d+$", "", rownames(neg_measurement))
pos_measurement$Group <- gsub("\\d+$", "", rownames(pos_measurement))
mix_measurement$Group <- gsub("\\d+$", "", rownames(mix_measurement))

# PCA
neg_pca <- prcomp(neg_measurement[, 1:ncol(neg_measurement)-1], scale = TRUE)
pos_pca <- prcomp(pos_measurement[, 1:ncol(pos_measurement)-1], scale = TRUE)
mix_pca <- prcomp(mix_measurement[, -18318], scale = TRUE)

# # plot the PCA
# fviz_pca_ind(neg_pca, 
#              geom = "point", 
#              palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00AFBB", "#E7B800", "#FC4E07")) +
#   geom_mark_ellipse(aes(fill = neg_measurement$Group)) +
#   #coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
#   theme_minimal() +
#   ggtitle("PCA")
# 
# 
# fviz_pca_ind(mix_pca, 
#              geom = "point", 
#              palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00AFBB", "#E7B800", "#FC4E07")) +
#   geom_mark_ellipse(aes(fill = mix_measurement$Group)) +
#   #coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
#   theme_minimal() +
#   ggtitle("PCA")



# get the new dataframe with the mean value of each group
## Mean
neg_measurement_mean <- data.frame(
  Con. = rowMeans(neg_measurement[, 2:4]),
  InApop. = rowMeans(neg_measurement[, 5:7]),
  ExApop. = rowMeans(neg_measurement[, 8:10]),
  Fer. = rowMeans(neg_measurement[, 11:13]),
  Nec. = rowMeans(neg_measurement[, 14:16]),
  row.names = neg_measurement$feature
  )
pos_measurement_mean <- data.frame(
  Con. = rowMeans(pos_measurement[, 2:4]),
  InApop. = rowMeans(pos_measurement[, 5:7]),
  ExApop. = rowMeans(pos_measurement[, 8:10]),
  Fer. = rowMeans(pos_measurement[, 11:13]),
  Nec. = rowMeans(pos_measurement[, 14:16]),
  row.names = pos_measurement$feature
  )
## all
neg_measurement_all <- data.frame(
  Con. = neg_measurement[, 2:4],
  InApop. = neg_measurement[, 5:7],
  ExApop. = neg_measurement[, 8:10],
  Fer. = neg_measurement[, 11:13],
  Nec. = neg_measurement[, 14:16], 
  row.names = neg_measurement$feature
  )
colnames(neg_measurement_all) <- c("Con.1", "Con.2", "Con.3", 
                                   "InApop.1", "InApop.2", "InApop.3", 
                                   "ExApop.1", "ExApop.2", "ExApop.3", 
                                   "Fer.1", "Fer.2", "Fer.3", 
                                   "Nec.1", "Nec.2", "Nec.3")
pos_measurement_all <- data.frame(
  Con. = pos_measurement[, 2:4],
  InApop. = pos_measurement[, 5:7],
  ExApop. = pos_measurement[, 8:10],
  Fer. = pos_measurement[, 11:13],
  Nec. = pos_measurement[, 14:16], 
  row.names = pos_measurement$feature
  )
colnames(pos_measurement_all) <- c("Con.1", "Con.2", "Con.3", 
                                   "InApop.1", "InApop.2", "InApop.3", 
                                   "ExApop.1", "ExApop.2", "ExApop.3", 
                                   "Fer.1", "Fer.2", "Fer.3", 
                                   "Nec.1", "Nec.2", "Nec.3")

# get the mix dataframe with positive and negative mode
mix_measurement_mean <- rbind(neg_measurement_mean, pos_measurement_mean)
mix_measurement_all <- rbind(neg_measurement_all, pos_measurement_all)

#  Z-score calculation
## mean
neg_measurement_z <- t(scale(t(neg_measurement_mean)))
pos_measurement_z <- t(scale(t(pos_measurement_mean)))
mix_measurement_z <- t(scale(t(mix_measurement_mean)))
## all
neg_measurement_all_z <- t(scale(t(neg_measurement_all)))
pos_measurement_all_z <- t(scale(t(pos_measurement_all)))
mix_measurement_all_z <- t(scale(t(mix_measurement_all)))


# correlation calculation + plotting
neg_mean_cor <- cor(neg_measurement_z)
pos_mean_cor <- cor(pos_measurement_z)
mix_mean_cor <- cor(mix_measurement_z)

col <- colorRampPalette(c("navy", "#FFFFFF", "#ee3430"))(200)

## mean
png("Ploting/1_output/NEG_Correlation.png", width = 1600, height = 1200, res = 300)
corrplot(neg_mean_cor, method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

png("Ploting/1_output/POS_Correlation.png", width = 1600, height = 1200, res = 300)
corrplot(pos_mean_cor, method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

png("Ploting/1_output/Mix_Correlation.png", width = 1600, height = 1200, res = 300)
corrplot(mix_mean_cor, method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

## all
png("Ploting/1_output/NEG_Correlation_all.png", width = 1800, height = 1400, res = 300)
corrplot(cor(neg_measurement_all_z), method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

png("Ploting/1_output/POS_Correlation_all.png", width = 1800, height = 1400, res = 300)
corrplot(cor(pos_measurement_all_z), method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

png("Ploting/1_output/Mix_Correlation_all.png", width = 1800, height = 1400, res = 300)
corrplot(cor(mix_measurement_all_z), method = "color", tl.col = "black", col = col, tl.srt = 45)
dev.off()

# PCA
## calculation
### mean
neg_measurement_pca <- prcomp(neg_measurement_z)
pos_measurement_pca <- prcomp(pos_measurement_z)
mix_measurement_pca <- prcomp(mix_measurement_z)

### all
neg_measurement_all_pca <- prcomp(neg_measurement_all_z)
pos_measurement_all_pca <- prcomp(pos_measurement_all_z)
mix_measurement_all_pca <- prcomp(mix_measurement_all_z)

## Scree plot
### mean
factoextra::fviz_eig(neg_measurement_pca)
factoextra::fviz_eig(pos_measurement_pca)
factoextra::fviz_eig(mix_measurement_pca)

### all
factoextra::fviz_eig(neg_measurement_all_pca)
factoextra::fviz_eig(pos_measurement_all_pca)
factoextra::fviz_eig(mix_measurement_all_pca)

## PCA plot
### extract the PCA matrix
#### mean
neg_mean_pca_rotation_wide <- neg_measurement_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")
pos_mean_pca_rotation_wide <- pos_measurement_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")
mix_mean_pca_rotation_wide <- mix_measurement_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")

#### all
neg_all_pca_rotation_wide <- neg_measurement_all_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")
neg_all_pca_rotation_wide$Group <- gsub("\\d+$", "", neg_all_pca_rotation_wide$Group)
neg_all_pca_rotation_wide$Group <- factor(neg_all_pca_rotation_wide$Group, levels = c("Con.", "InApop.", "ExApop.", "Fer.", "Nec."))
pos_all_pca_rotation_wide <- pos_measurement_all_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")
pos_all_pca_rotation_wide$Group <- gsub("\\d+$", "", pos_all_pca_rotation_wide$Group)
pos_all_pca_rotation_wide$Group <- factor(pos_all_pca_rotation_wide$Group, levels = c("Con.", "InApop.", "ExApop.", "Fer.", "Nec."))
mix_all_pca_rotation_wide <- mix_measurement_all_pca$rotation %>% 
  as.data.frame() %>%
  rownames_to_column("Group")
mix_all_pca_rotation_wide$Group <- gsub("\\d+$", "", mix_all_pca_rotation_wide$Group)
mix_all_pca_rotation_wide$Group <- factor(mix_all_pca_rotation_wide$Group, levels = c("Con.", "InApop.", "ExApop.", "Fer.", "Nec."))

### plot
#### mean
ggplot(mix_mean_pca_rotation_wide, aes(x = PC1, y = PC2, color = Group))+
  geom_point()+
  theme_classic()

#### all
ggbiplot(mix_measurement_all_pca,
         ) + 
  labs(fill = "Group", color = "Group")






ggplot(mix_all_pca_rotation_wide, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.5) +
  geom_encircle(aes(x = PC1, y = PC2, color = Group), 
                data = mix_all_pca_rotation_wide, 
                size = 2, 
                expand = 0.05) +
  coord_fixed(xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7)) +
  labs(x = "PC1 (46.9%)", y = "PC2 (12.4%)") +
  theme_classic() +
  theme(
    legend.key.size = unit(0.5, "cm"),
    legend.position = c(0.1, 0.87),
    legend.title = element_blank()
  )

fviz_pca_var(mix_measurement_all_pca)



