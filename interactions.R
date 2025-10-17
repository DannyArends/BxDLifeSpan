#
# lifespan.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

library(lme4)
library(parallel)

source("shared.R")

### individual Level mapping (Random effects used to deal with correlation structure caused by multiple strain observations)

n.cores <- (detectCores() / 2) - 1
clust <- makeCluster(n.cores)
clusterExport(clust, "ldata")
clusterExport(clust, "geno")
clusterCall(clust, function() library(lme4))

pvalsL <- parLapply(clust, seq(30, 780, 30), function(day){
  isDIET <- ldata[which(grepl("BXD", ldata[, "StrainName"]) & ldata[, "AgeAtDeath..days."] >= day),]
  isDIET <- isDIET[which(isDIET[, "StrainName"] %in% colnames(geno)),]
  genoS <- geno[, isDIET[, "StrainName"]]
  strains <- isDIET[, "StrainName"]

  Y <- isDIET[, "AgeAtDeath..days."]
  AgeAtSetUp <- isDIET[, "AgeAtSetUp.in.colony..days."]
  vivarium <- as.factor(isDIET[, "vivarium"])
  diet <- as.factor(isDIET[, "Diet"])

  null <- lmer(Y ~ diet + vivarium + AgeAtSetUp + (1 | strains), REML = FALSE)

  pvals.main <- c()
  pvals.int <- c()
  for(x in c(1:nrow(geno))) {
    gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2

    fd <- data.frame(Y, strains, diet, AgeAtSetUp, vivarium, gts)
    full <- lmer(Y ~ diet + vivarium + AgeAtSetUp + gts + (1 | strains), data = fd, REML = FALSE)
    full.int <- lmer(Y ~ diet + vivarium + AgeAtSetUp + gts + diet:gts + (1 | strains), data = fd, REML = FALSE)

    iix <- which(apply(apply(fd, 1, is.na),2,sum) > 0)

    if(length(iix) > 0) {
      ## If we have missing genotypes, we need to account for this, specify an ALT model
      fd <- fd[-iix,]
      alt0 <- lmer(Y ~ diet + vivarium + AgeAtSetUp + (1 | strains), data = fd, REML = FALSE)
      pvals.main <- c(pvals.main, as.numeric(na.omit(anova(alt0, full)[, "Pr(>Chisq)"])))
      pvals.int <- c(pvals.int, as.numeric(na.omit(anova(full, full.int)[, "Pr(>Chisq)"])))
    }else{
      ## No missing genotypes, so just use the global NULL model which has all data
      pvals.main <- c(pvals.main, as.numeric(na.omit(anova(null, full)[, "Pr(>Chisq)"])))
      pvals.int <- c(pvals.int, as.numeric(na.omit(anova(full, full.int)[, "Pr(>Chisq)"])))
    }
    if(x %% 100 == 1) cat(x, "\n")
  }
  return(cbind(main = pvals.main, int = pvals.int))
})
stopCluster(clust)

pvalM.main <- c()
pvalM.int <- c()
for(x in 1:length(pvalsL)){ 
  pvalM.main <- cbind(pvalM.main, pvalsL[[x]][,1])
  pvalM.int <- cbind(pvalM.int, pvalsL[[x]][,2])
}

res.main <- cbind(map, -log10(pvalM.main))
colnames(res.main) <- c("Chr","Locus","cM","Mb", seq(30, 780, 30))

res.int <- cbind(map, -log10(pvalM.int))
colnames(res.int) <- c("Chr","Locus","cM","Mb", seq(30, 780, 30))

write.table(res.main, file = paste0("output/main_ind_Progressive.txt"), sep = "\t", quote = FALSE)
write.table(res.int, file = paste0("output/int_ind_Progressive.txt"), sep = "\t", quote = FALSE)

op <- par(mfrow = c(2,1))
image(-log10(pvalM.main), breaks = c(0,1,2.5,100), col = c("white", "#E1E1E1", "red"), xaxt="n", yaxt = "n")
plot(-log10(pvalM.main[,1]), col = as.numeric(as.factor(map[, "Chr"])), xaxs="i")

op <- par(mfrow = c(2,1))
image(-log10(pvalM.int), breaks = c(0,1,2.5,100), col = c("white", "#E1E1E1", "red"), xaxt="n", yaxt = "n")
plot(-log10(pvalM.int[,1]), col = as.numeric(as.factor(map[, "Chr"])), xaxs="i")

image(-log10(pvalM.int[which(map[, "Chr"] == 1),]), breaks = c(0,1,2.5,100), col = c("white", "#E1E1E1", "red"), xaxt="n", yaxt = "n")



