#
# Lifespan - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
setwd("C:/Github/BxDLifespan/")
ldata <- read.table("data/BxD_Lifespan.txt", sep = "\t", header = TRUE, row.names=2)
bdata <- read.table("data/BxD_Bodyweights.txt", sep = "\t", header = TRUE)

### One eartag is found attached to two different mice
eartag <- names(which(table(bdata[,1]) > 1))
bdata <- bdata[-which(bdata[,1] == eartag),]
ldata <- ldata[-which(rownames(ldata) == eartag),]

### Only take animals which are not "Harvested"
iix <- which(bdata[, "Status"] %in% c("Harvest", "Other"))
bdata <- bdata[-iix,]

### Still some missing mice (Yes bw, no lifespan)
rownames(bdata) <- bdata[,1]

### Mice which are matched
bdata <- bdata[which(rownames(bdata) %in% rownames(ldata)),]
ldata <- ldata[which(rownames(ldata) %in% rownames(bdata)),]

### Sort based on the order in lifespan

bdata <- bdata[rownames(ldata),]

### Load in genotypes
geno <- read.table("http://files.genenetwork.org/current/GN600/BXD_current_rev050423.geno", 
                   sep = "\t", skip=23, header = TRUE, na.strings = c("U", "", "NA"))
map <- geno[, 1:4]
geno <- geno[, -c(1:4)]
rownames(geno) <- map[,2]
rownames(map) <- map[,2]

op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF")){
  isCD <- ldata[which(ldata[, "Diet"] == diet),]
  mtable <- table(isCD[, "StrainName"])
  mtable <- mtable[which(mtable >= 3)]
  strains <- names(mtable)
  bxds <- strains[grep("BXD", strains)]

  strainM <- round(unlist(lapply(bxds, function(x){ mean(isCD[which(isCD[, "StrainName"] == x), "AgeAtDeath..days."]); })), 0)
  names(strainM) <- bxds

  ### Mapping Strain means
  genoS <- geno[,names(strainM)]

  pvals <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
    pvals <- c(pvals, anova(lm(strainM ~ gts))[[5]][1])
    if(x %% 1000 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0(diet, "_strainMeans.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

library(lme4)
### individual Level mapping (Random effect)
op <- par(mfrow = c(2,1))
for(diet in c("CD", "HF")){
  isCD <- ldata[which(ldata[, "Diet"] == diet & grepl("BXD", ldata[, "StrainName"])),]
  isCD <- isCD[which(isCD[, "StrainName"] %in% colnames(geno)),]
  genoS <- geno[, isCD[, "StrainName"]]
  strains <- isCD[, "StrainName"]
  null <- lmer(isCD[, "AgeAtDeath..days."] ~ 1 + (1 | strains), REML = FALSE)

  pvals <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
    full <- lmer(isCD[, "AgeAtDeath..days."] ~ gts + (1 | strains), REML = FALSE)
    if(any(is.na(gts))){
      ## If we have missing genotypes, we need to account for this, specify an ALT model
      alt <- lmer(isCD[which(!is.na(gts)), "AgeAtDeath..days."] ~ 1 + (1 | strains[which(!is.na(gts))]), REML = FALSE)
      pvals <- c(pvals, as.numeric(na.omit(anova(alt, full)[, "Pr(>Chisq)"])))
    }else{
      ## No missing genotypes, so just use the global NULL model which has all data
      pvals <- c(pvals, as.numeric(na.omit(anova(null, full)[, "Pr(>Chisq)"])))
    }
    if(x %% 100 == 1) cat(x, "\n")
  }
  write.table(cbind(map, -log10(pvals)), file = paste0(diet, "_indWeighted.txt"), sep = "\t", quote = FALSE)
  plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])))
}

### CTL mapping

### Interactions together next Wednesday



