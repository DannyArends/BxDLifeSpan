#
# bodyweight.R Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
# AgeAtSetUp == DateWeightBaseline == DateDietStart 
# Weight00.Baseline ~~ AgeAtSetUp


library(lme4)
library(parallel)

source("shared.R")

### Mapping Weighted Strain means
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF")) {
  isDIET <- bdata[which(ldata[, "Diet"] == diet),]
  mtable <- table(isDIET[, "StrainName"])
  mtable <- mtable[which(mtable >= 3)]
  strains <- names(mtable)
  bxds <- strains[grep("BXD", strains)]

  weights <- mtable[bxds]

  strainM <- round(unlist(lapply(bxds, function(x){ 
    mean(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtDeath..days."]); })), 2)
  names(strainM) <- bxds

  genoS <- geno[,names(strainM)]

  pvals <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
    pvals <- c(pvals, anova(lm(strainM ~ gts, weights = weights))[[5]][1])
    if(x %% 1000 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0("output/", diet, "_wMeans.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

