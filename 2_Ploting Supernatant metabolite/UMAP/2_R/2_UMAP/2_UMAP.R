#----import data----
mix_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/Metaboanalyst/Mix/Mix_iden.csv", stringsAsFactors = FALSE)

POS_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_POS/3_QI Report/POS_Compound Measurement.csv", stringsAsFactors = FALSE)

NEG_measure <- read.csv("1_Data/20231204 Cell-Death Metabolite_Sup/20231204_NEG/3_QI Report/NEG_Compound Measurement.csv", stringsAsFactors = FALSE)

#----library----
library(tidyverse)
library(umap)

#----data processing----
mix_measure_t <- t(mix_measure)

pos_measure <- POS_measure[, 14:28]
pos_measure_t <-  t(pos_measure)

neg_measure <- NEG_measure[, 14:28]
neg_measure_t <- t(neg_measure)

##----rename the colnames----
colnames(mix_measure_t) <- mix_measure_t[1,]
colnames(pos_measure_t) <- POS_measure[, 1]
colnames(neg_measure_t) <- NEG_measure[, 1]

mix_measure_t <- mix_measure_t[-1,] # remove the first row

##----convert to dataframe----
mix_measure_t <- data.frame(mix_measure_t)
pos_measure_t <- data.frame(pos_measure_t)
neg_measure_t <- data.frame(neg_measure_t)

##----convert to numeric, for umap data----
mix_measure_data <- mutate_all(mix_measure_t, function(x) as.numeric(as.character(x)))
pos_measure_data <- mutate_all(pos_measure_t, function(x) as.numeric(as.character(x)))
neg_measure_data <- mutate_all(neg_measure_t, function(x) as.numeric(as.character(x)))

##----convert to factor, for umap label----
mix_measure_label <- as.factor(rownames(mix_measure_t))
pos_measure_label <- as.factor(rownames(mix_measure_t))
neg_measure_label <- as.factor(rownames(mix_measure_t))

#----umap function----
plot.umap <- function(x, labels,
                      main="A UMAP visualization of the Cell-Death Metabolite Dataset (POS +NEG)",
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


#----umap calculation----
mix_umap <- umap(mix_measure_data)
pos_umap <- umap(pos_measure_data)
neg_umap <- umap(neg_measure_data)

#----plot----
set.seed(1)
mix.umap.plot <- plot.umap(mix_umap, mix_measure_label)
plot.map(pos_umap, pos_measure_label)
plot.umap(neg_umap, neg_measure_label)

#----save----
dev.print(png, file = "1_Data/20231204 Cell-Death Metabolite_Sup/umap/Mix_umap.png", )