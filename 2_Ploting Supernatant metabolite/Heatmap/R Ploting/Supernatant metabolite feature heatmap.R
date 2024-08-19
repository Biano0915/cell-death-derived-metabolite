library(tidyverse)
library(readxl)
library(pheatmap)

# Load the data
neg_feature <- read_excel("R Ploting/Feature.xlsx", sheet = "Negative", col_names = TRUE)
neg_feature_to_compound_1 <- read_excel("R Ploting/Feature.xlsx", sheet = "neg_feature_to_compound_1", col_names = TRUE)
neg_feature_to_compound_2 <- read_excel("R Ploting/Feature.xlsx", sheet = "neg_feature_to_compound_2", col_names = TRUE)

pos_feature <- read_excel("R Ploting/Feature.xlsx", sheet = "Positive", col_names = TRUE)
pos_feature_to_compound_1 <- read_excel("R Ploting/Feature.xlsx", sheet = "pos_feature_to_compound_1", col_names = TRUE)
pos_feature_to_compound_2 <- read_excel("R Ploting/Feature.xlsx", sheet = "pos_feature_to_compound_2", col_names = TRUE)

# rename the columns
colnames(neg_feature) <- c("feature", "Con.1", "Con.2", "Con.3", 
                           "InApop.1", "InApop.2", "InApop.3",
                           "ExApop.1", "ExApop.2", "ExApop.3",
                           "Fer.1", "Fer.2", "Fer.3",
                           "Nec.1", "Nec.2", "Nec.3")
colnames(pos_feature) <- c("feature", "Con.1", "Con.2", "Con.3", 
                           "InApop.1", "InApop.2", "InApop.3",
                           "ExApop.1", "ExApop.2", "ExApop.3",
                           "Fer.1", "Fer.2", "Fer.3",
                           "Nec.1", "Nec.2", "Nec.3")

# convert to a frame frame and set row names
neg_feature <- as.data.frame(neg_feature)
row.names(neg_feature) <- neg_feature$feature
neg_feature <- neg_feature[, -1]

pos_feature <- as.data.frame(pos_feature)
row.names(pos_feature) <- pos_feature$feature
pos_feature <- pos_feature[, -1]

# calculate the mean of each treatment
neg_feature_mean <- data.frame(
  Con. = rowMeans(neg_feature[, 1:3]),
  InApop. = rowMeans(neg_feature[, 4:6]),
  ExApop. = rowMeans(neg_feature[, 7:9]),
  Fer. = rowMeans(neg_feature[, 10:12]),
  Nec. = rowMeans(neg_feature[, 13:15])
)

pos_feature_mean <- data.frame(
  Con. = rowMeans(pos_feature[, 1:3]),
  InApop. = rowMeans(pos_feature[, 4:6]),
  ExApop. = rowMeans(pos_feature[, 7:9]),
  Fer. = rowMeans(pos_feature[, 10:12]),
  Nec. = rowMeans(pos_feature[, 13:15])
)

# calculate the z-score
neg_feature_z <- t(scale(t(neg_feature)))
pos_feature_z <- t(scale(t(pos_feature)))

neg_feature_mean_z <- t(scale(t(neg_feature_mean)))
pos_feature_mean_z <- t(scale(t(pos_feature_mean)))


# setup the named vector
neg_compound_names <- setNames(neg_feature_to_compound_2$Compound, neg_feature_to_compound_2$Feature)
pos_compound_names <- setNames(pos_feature_to_compound_2$Compound, pos_feature_to_compound_2$Feature)

# rename the row names
rownames(neg_feature_z) <- neg_compound_names[rownames(neg_feature_z)]
rownames(pos_feature_z) <- pos_compound_names[rownames(pos_feature_z)]

rownames(neg_feature_mean_z) <- neg_compound_names[rownames(neg_feature_mean_z)]
rownames(pos_feature_mean_z) <- pos_compound_names[rownames(pos_feature_mean_z)]

# remove the unamed rows
neg_feature_z <- neg_feature_z[!is.na(rownames(neg_feature_z)), ]
pos_feature_z <- pos_feature_z[!is.na(rownames(pos_feature_z)), ]

neg_feature_mean_z <- neg_feature_mean_z[!is.na(rownames(neg_feature_z)), ]
pos_feature_mean_z <- pos_feature_mean_z[!is.na(rownames(pos_feature_z)), ]

# setup for heatmap
## calculate the data range
neg_feature_range <- range(neg_feature_z)
pos_feature_range <- range(pos_feature_z)

neg_feature_mean_range <- range(neg_feature_mean_z)
pos_feature_mean_range <- range(pos_feature_mean_z)

## create the symmetric breaks around zero
max_val_neg <- max(abs(neg_feature_range))
breaks_neg <- seq(-max_val_neg, max_val_neg, length.out = 100)
max_val_pos <- max(abs(pos_feature_range))
breaks_pos <- seq(-max_val_pos, max_val_pos, length.out = 100)

max_val_neg_mean <- max(abs(neg_feature_mean_range))
breaks_neg_mean <- seq(-max_val_neg_mean, max_val_neg_mean, length.out = 100)
max_val_pos_mean <- max(abs(pos_feature_mean_range))
breaks_pos_mean <- seq(-max_val_pos_mean, max_val_pos_mean, length.out = 100)

# plot the heatmap
pheatmap(neg_feature_z, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = breaks_neg,
         main = "Anionic Features", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE)
pheatmap(pos_feature_z,
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         breaks = breaks_pos,
         main = "Cationic Features", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE)


pheatmap(neg_feature_mean_z, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Anionic Features", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE)
pheatmap(pos_feature_mean_z, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Cationic Features", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE)

# save the z-score data
write.csv(neg_feature_z, "1_output/1_feature z score/neg_feature_z.csv")
write.csv(pos_feature_z, "1_output/1_feature z score/pos_feature_z.csv")


