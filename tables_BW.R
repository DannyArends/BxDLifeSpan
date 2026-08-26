#
# tables.R - Provide nice tables for top, flanking, and effect sizes for Vta & Soma Loci in BxD
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

source("shared.R")
wMeans.CD <- read.table("output/BW_CD_wMeans_Progressive.txt", sep = "\t")
wMeans.HF <- read.table("output/BW_HF_wMeans_Progressive.txt", sep = "\t")

mA3 <- sort(unique(as.numeric(unlist(apply(wMeans.CD[, -c(1:4)], 2, function(x){ which(x > 3) })))))

tpMax <- apply(wMeans.CD[mA3, -c(1:4)],1, which.max)

cat("", file = "wMean.summary.CD.txt")
for(i in 1:length(mA3)){
  chr <- wMeans.CD[mA3[i], 1]
  locus <- wMeans.CD[mA3[i], 2]
  mb <- wMeans.CD[mA3[i], 4]
  tp <- colnames(wMeans.CD)[-c(1:4)][tpMax[i]]
  lod <- wMeans.CD[mA3[i], -c(1:4)][tpMax[i]]
  cat(locus, "\t", chr, "\t", mb, "\t", gsub("X", "", tp), "\t", as.numeric(lod), "\n", file = "wMean.BW.summary.CD.txt", append = TRUE)
}

