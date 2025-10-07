#
# lifespan.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow

library(lme4)
library(parallel)

source("shared.R")

### Mapping Weighted Strain means
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF")) {
  isDIET <- ldata[which(ldata[, "Diet"] == diet),]
  mtable <- table(isDIET[, "StrainName"])
  mtable <- mtable[which(mtable >= 3)]
  strains <- names(mtable)
  bxds <- strains[grep("BXD", strains)]

  weights <- mtable[bxds]

  strainM <- round(unlist(lapply(bxds, function(x){ mean(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtDeath..days."]); })), 0)
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

### individual Level mapping (Random effects used to deal with correlation structure caused by multiple strain observations)
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF")){
  isDIET <- ldata[which(ldata[, "Diet"] == diet & grepl("BXD", ldata[, "StrainName"])),]
  isDIET <- isDIET[which(isDIET[, "StrainName"] %in% colnames(geno)),]
  genoS <- geno[, isDIET[, "StrainName"]]
  strains <- isDIET[, "StrainName"]
  null <- lmer(isDIET[, "AgeAtDeath..days."] ~ 1 + (1 | strains), REML = FALSE)

  pvals <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
    full <- lmer(isDIET[, "AgeAtDeath..days."] ~ gts + (1 | strains), REML = FALSE)
    if(any(is.na(gts))){
      ## If we have missing genotypes, we need to account for this, specify an ALT model
      alt <- lmer(isDIET[which(!is.na(gts)), "AgeAtDeath..days."] ~ 1 + (1 | strains[which(!is.na(gts))]), REML = FALSE)
      pvals <- c(pvals, as.numeric(na.omit(anova(alt, full)[, "Pr(>Chisq)"])))
    }else{
      ## No missing genotypes, so just use the global NULL model which has all data
      pvals <- c(pvals, as.numeric(na.omit(anova(null, full)[, "Pr(>Chisq)"])))
    }
    if(x %% 100 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0("output/", diet, "_ind.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

### Progressive Mapping, Weighted Strain means
for(diet in c("CD", "HF")) {
  pvalsM <- c()

  n.cores <- (detectCores() / 2) - 1
  clust <- makeCluster(n.cores)
  clusterExport(clust, "ldata")
  clusterExport(clust, "diet")
  clusterExport(clust, "geno")

  pvalsL <- parLapply(clust, seq(30, 780, 30), function(day){
    isDIET <- ldata[which(ldata[, "Diet"] == diet & ldata[, "AgeAtDeath..days."] >= day),]
    mtable <- table(isDIET[, "StrainName"])
    mtable <- mtable[which(mtable >= 3)]
    strains <- names(mtable)
    bxds <- strains[grep("BXD", strains)]
    weights <- mtable[bxds]

    strainM <- round(unlist(lapply(bxds, function(x){ mean(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtDeath..days."]); })), 0)
    names(strainM) <- bxds

    genoS <- geno[,names(strainM)]

    pvals <- c()
    for(x in c(1:nrow(geno))) {
      gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
      pvals <- c(pvals, anova(lm(strainM ~ gts, weights = weights))[[5]][1])
    }
    return(pvals)
  })
  stopCluster(clust)
  plot(-log10(pvalsL[[1]]), col = as.numeric(as.factor(map[, "Chr"])))

  pvalM <- c()
  for(x in 1:length(pvalsL)){ pvalM <- cbind(pvalM, pvalsL[[x]]) }

  op <- par(mfrow=c(2,1))
  image(-log10(pvalM), breaks = c(0,1,2,3,4,5,6,10), 
    col = c("white", "#E1E1E1", "#D0D0D0", "#00C0C0","#00B0B0","#00A0A0","#009090"), xaxt="n", yaxt = "n")
  axis(2, at = (1:15 - 0.5) / 15, seq(30, 900, 60), las = 2)
  plot(-log10(pvalsL[[1]]), col = as.numeric(as.factor(map[, "Chr"])), xaxs="i")

  write.table(cbind(map, -log10(pvalM)), file = paste0("output/", diet, "_wMeans_Progressive.txt"), sep = "\t", quote = FALSE)
}

### Interactions together next Wednesday



