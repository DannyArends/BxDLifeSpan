
library(lme4)
library(parallel)

source("shared.R")

# Cut-AgeAtSetup
ldata[which(ldata[,"AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
ldata <- cbind(ldata, AgeAtSetUpGroup = cut(as.numeric(ldata[, "AgeAtSetUp.in.colony..days."]), seq(0, 1000, 180)))
ldata <- ldata[-which(is.na(ldata[, "AgeAtSetUpGroup"])),]

### Mapping Weighted Strain means (interactions)
op <- par(mfrow = c(2,1))
diets <- c("HF", "NAM")

isDIET <- ldata[which(ldata[, "Diet"] %in% diets & grepl("BXD", ldata[, "StrainName"])),]
isDIET <- isDIET[which(isDIET[, "StrainName"] %in% colnames(geno)),]
genoS <- geno[, isDIET[, "StrainName"]]
strains <- isDIET[, "StrainName"]

Y <- isDIET[, "AgeAtDeath..days."]
AgeAtSetUp <- as.numeric(isDIET[, "AgeAtSetUp.in.colony..days."])
AgeAtSetUp <- as.factor(cut(AgeAtSetUp, seq(0, 1000, 180))) ### Threshold for less extreme

diet <- isDIET[, "Diet"]

pvals <- c()
for(x in 1:nrow(geno)) {
  gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2

  fd <- cbind(Y, AgeAtSetUp, diet, gts)

  null <- lmer(Y ~ AgeAtSetUp + diet + gts + (1 | strains), REML = FALSE)
  full <- lmer(Y ~ AgeAtSetUp + diet + gts + gts:diet + (1 | strains), REML = FALSE)

  iix <- which(apply(apply(fd, 1, is.na),2,sum) > 0)

  if(length(iix) > 0) {
    ## If we have missing genotypes, we need to account for this, specify an ALT model
    alt0 <- lmer(Y[-iix] ~ AgeAtSetUp[-iix] + diet[-iix] + gts[-iix] + (1 | strains[-iix]), REML = FALSE)
    pvals <- c(pvals, as.numeric(na.omit(anova(alt0, full)[, "Pr(>Chisq)"])))
  }else{
    ## No missing genotypes, so just use the global NULL model which has all data
    pvals <- c(pvals, as.numeric(na.omit(anova(null, full)[, "Pr(>Chisq)"])))
  }
  if(x %% 100 == 1) cat(x, "\n")
}
write.table(cbind(map, -log10(pvals)), file = paste0("output/NAM_Interaction.txt"), sep = "\t", quote = FALSE)


edge1 <- "rsm10000000387"
edge2 <- "rsm10000000387"
