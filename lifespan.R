#
# lifespan.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

library(lme4)
library(parallel)

source("shared.R")

# Cut-AgeAtSetup
ldata[which(ldata[,"AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
ldata <- cbind(ldata, AgeAtSetUpGroup = cut(as.numeric(ldata[, "AgeAtSetUp.in.colony..days."]), seq(0, 1000, 180)))
ldata <- ldata[-which(is.na(ldata[, "AgeAtSetUpGroup"])),]
write.table(ldata, file = "output/BxDLifespan_SetUpAge.txt", sep = "\t", quote = FALSE)

### Write out the animals we used:
write.table(cbind(rownames(ldata), ldata[, "Diet"], ldata[, "StrainName"]), file = "output/EarTagNumberCurrent.txt", sep = "\t", quote = FALSE, row.names=FALSE)

### Mapping Weighted Strain means
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF", "NAM")) {
  isDIET <- ldata[which(ldata[, "Diet"] == diet),]
  mtable <- table(isDIET[, "StrainName"])
  mtable <- mtable[which(mtable >= 3)]
  strains <- names(mtable)
  bxds <- strains[grep("BXD", strains)]

  weights <- mtable[bxds]

  strainM <- round(unlist(lapply(bxds, function(x){ mean(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtDeath..days."]); })), 0)
  names(strainM) <- bxds

  ## AgeAtSetUp P ~ 0.1
  AgeAtSetUp <- round(unlist(lapply(bxds, function(x){ mean(as.numeric(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtSetUp.in.colony..days."])); })), 0)
  AgeAtSetUp <- as.factor(cut(AgeAtSetUp, seq(0, 1000, 180))) ### Threshold for less extreme
  names(AgeAtSetUp) <- bxds

  genoS <- geno[,names(strainM)]

  pvals <- c()
  for(x in c(1:nrow(geno))) {
      gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
      pvals <- c(pvals, anova(lm(strainM ~ AgeAtSetUp + gts, weights = weights))[[5]][2])
    if(x %% 1000 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0("output/", diet, "_wMeans.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

### individual Level mapping (Random effects used to deal with correlation structure caused by multiple strain observations)
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF", "NAM")){
  isDIET <- ldata[which(ldata[, "Diet"] == diet & grepl("BXD", ldata[, "StrainName"])),]
  isDIET <- isDIET[which(isDIET[, "StrainName"] %in% colnames(geno)),]
  genoS <- geno[, isDIET[, "StrainName"]]
  strains <- isDIET[, "StrainName"]

  Y <- isDIET[, "AgeAtDeath..days."]
  AgeAtSetUp <- as.numeric(isDIET[, "AgeAtSetUp.in.colony..days."])
  AgeAtSetUp <- as.factor(cut(AgeAtSetUp, seq(0, 1000, 180))) ### Threshold for less extreme
  vivarium <- as.factor(isDIET[, "vivarium"])

  aa <- boxplot(Y ~ AgeAtSetUp)
  pdf(paste0("output/boxplot_AgeAtSetup_",diet,".pdf"), width = 24, height = 18)
  plot(Y ~ AgeAtSetUp)
  text(1:length(aa$n), rep(aa$stats[3,]-50, length(aa$n)), paste0("n=", aa$n), col = "white")
  dev.off()

  null <- lmer(Y ~ AgeAtSetUp + (1 | strains), REML = FALSE)
  if(length(levels(vivarium)) >1) null <- lmer(Y ~ vivarium + AgeAtSetUp + (1 | strains), REML = FALSE)

  pvals <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2

    fd <- cbind(Y, AgeAtSetUp, vivarium, gts)

    if(length(levels(vivarium)) >1){
      full <- lmer(Y ~ vivarium + AgeAtSetUp + gts + (1 | strains), REML = FALSE)
    }else{
      full <- lmer(Y ~ AgeAtSetUp + gts + (1 | strains), REML = FALSE)
    }

    iix <- which(apply(apply(fd, 1, is.na),2,sum) > 0)

    if(length(iix) > 0) {
      ## If we have missing genotypes, we need to account for this, specify an ALT model
     if(length(levels(vivarium)) >1){
        alt0 <- lmer(Y[-iix] ~ vivarium[-iix] + AgeAtSetUp[-iix] + (1 | strains[-iix]), REML = FALSE)
      }else{
        alt0 <- lmer(Y[-iix] ~ AgeAtSetUp[-iix] + (1 | strains[-iix]), REML = FALSE)
      }
      #alt0 <- lmer(Y[-iix] ~ (1 | strains[-iix]), REML = FALSE)
      pvals <- c(pvals, as.numeric(na.omit(anova(alt0, full)[, "Pr(>Chisq)"])))
    }else{
      ## No missing genotypes, so just use the global NULL model which has all data
      pvals <- c(pvals, as.numeric(na.omit(anova(null, full)[, "Pr(>Chisq)"])))
    }
    if(x %% 100 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0("output/", diet, "_ind.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

pvalsILM <- pvals
#### UPDATe for NAME mapping, remove the vivarium from the model , when mapping NAMs

### Progressive Mapping, Weighted Strain means
for(diet in c("CD", "HF")) {
  pvalsM <- c()
  nVector <- c()

  n.cores <- (detectCores() / 2) - 1
  clust <- makeCluster(n.cores)
  clusterExport(clust, "ldata")
  clusterExport(clust, "diet")
  clusterExport(clust, "geno")

  pvalsL <- parLapply(clust, seq(20, 780, 30), function(day){
    isDIET <- ldata[which(ldata[, "Diet"] == diet & ldata[, "AgeAtDeath..days."] >= day),]
    mtable <- table(isDIET[, "StrainName"])
    mtable <- mtable[which(mtable >= 3)]
    strains <- names(mtable)
    bxds <- strains[grep("BXD", strains)]
    weights <- mtable[bxds]

    strainM <- round(unlist(lapply(bxds, function(x){ mean(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtDeath..days."]); })), 0)
    names(strainM) <- bxds

    ## AgeAtSetUp P ~ 0.1
    AgeAtSetUp <- round(unlist(lapply(bxds, function(x){ mean(as.numeric(isDIET[which(isDIET[, "StrainName"] == x), "AgeAtSetUp.in.colony..days."])); })), 0)
    AgeAtSetUp <- as.factor(cut(AgeAtSetUp, seq(0, 1000, 180))) ### Threshold for less extreme
    names(AgeAtSetUp) <- bxds

    genoS <- geno[,names(strainM)]

    pvals <- c()
    for(x in c(1:nrow(geno))) {
      gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
      pvals <- c(pvals, anova(lm(strainM ~ AgeAtSetUp + gts, weights = weights))[[5]][2])
    }
    return(list(pvals, length(names(strainM))))
  })
  stopCluster(clust)
  ns <- unlist(lapply(pvalsL, "[", 2))
  names(ns) <- seq(20, 780, 30)

  pvalM <- c()
  for(x in 1:length(pvalsL)){ pvalM <- cbind(pvalM, pvalsL[[x]][[1]]) }

  op <- par(mfrow=c(2,1))
  image(-log10(pvalM), breaks = c(0,1,2,3,4,5,6,10), 
    col = c("white", "#E1E1E1", "#D0D0D0", "#00C0C0","#00B0B0","#00A0A0","#009090"), xaxt="n", yaxt = "n")
  #axis(2, at = (1:15 - 0.5) / 15, seq(30, 780, 30), las = 2)

  res <- cbind(map, -log10(pvalM))
  colnames(res) <- c("Chr","Locus","cM","Mb", seq(20, 780, 30))
  write.table(res, file = paste0("output/", diet, "_wMeans_Y_Progressive.txt"), sep = "\t", quote = FALSE)
  write.table(ns, file = paste0("output/", diet, "_wMeans_Y_Progressive_N.txt"), sep = "\t", quote = FALSE)
}

### Interactions together next Wednesday

### Progressive Mapping, Weighted Strain means

