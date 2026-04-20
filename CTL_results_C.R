#
# lifespan.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

library(RColorBrewer)
source("shared.R")

res.main <- -log10(read.table("output/CTL_CD_Pmatrix.txt", sep = "\t"))
res.int <- -log10(read.table("output/CTL_HF_Pmatrix.txt", sep = "\t"))

pdf(paste0("output/CTL_Results_Cor_at_TP.pdf"), width = 24, height = 11)
  op <- par(mfrow=c(1, 2))
  plot(x = as.numeric(gsub("X", "",colnames(cor(res.main[,"X470"], res.main))[-1])), 
       y=cor(res.main[,"X470"], res.main)[-1], main = "Correlation of CTL result TP to other TP (CD)", pch = 19, xlab = "Timepoint", ylab = "R")

  points(x = as.numeric(gsub("X", "",colnames(cor(res.main[,"X440"], res.main))[-1])), 
         y=cor(res.main[,"X440"], res.main)[-1], col = "green", pch = 19)
  legend("topleft", c("T440", "T470"), pch = 19, col = c("black", "green"))


  plot(x = as.numeric(gsub("X", "",colnames(cor(res.int[,"X470"], res.int))[-1])), 
       y=cor(res.int[,"X470"], res.int)[-1], main = "Correlation of CTL result TP to other TP (HF)", pch = 19, xlab = "Timepoint", ylab = "R")

  points(x = as.numeric(gsub("X", "",colnames(cor(res.int[,"X440"], res.int))[-1])), 
         y=cor(res.int[,"X440"], res.int)[-1], col = "green", pch = 19)

  legend("topleft", c("T440", "T470"), pch = 19, col = c("black", "green"))
dev.off()

