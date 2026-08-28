#
# effects.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

library(lme4)
library(RColorBrewer)
source("shared.R")

mdata <- read.table("output/table.txt", sep = "\t", header = TRUE)

lodfiles <- list(CD = "output/CD_wMeans_Y_Progressive.txt", HF = "output/HF_wMeans_Y_Progressive.txt")
lods <- lapply(lodfiles, read.table, sep = "\t", header = TRUE, check.names = FALSE)
lods <- lapply(lods, function(z){ rownames(z) <- z[, "Locus"]; z })
lodFloor <- -28

std <- function(x) sd(x)/sqrt(length(x))

for(x in 1:nrow(mdata)){
  chr <- mdata[x, "Chr"]
  top <- mdata[x, "Peak"]
  onchr <- map[map[, 1] == chr,]
  ii <- rownames(onchr)[which.min(abs(onchr[, 4] - top))]
  cat(mdata[x,1], chr, top, ii, "\n")
  pdf(paste0("output/effects/",mdata[x, 1],".pdf"), width = 12, height = 8)
  plot(c(20, 780), c(lodFloor, 20), t = "n", xlab = "Truncation Age [d]", ylab = "Actuarial Effect Size [d]", main = mdata[x, 1], xaxt="n", yaxt = "n")
  axis(1, at = seq(20, 780, 30), rep("", 26), tcl = -0.25)
  axis(1, at = seq(20, 780, 30)[seq(2, 26,2)], seq(20, 780, 30)[seq(2, 26,2)])
  axis(2, at = seq(-20, 20, 1), rep("", 41), tcl = -0.25)
  axis(2, at = seq(-20, 20, 5), seq(-20, 20, 5), las = 2)
  abline(h = -20, col = "gray")
  abline(h = lodFloor + 3, lty = 2, col = "gray")
  axis(4, at = lodFloor + seq(0, 8, 2), seq(0, 8, 2), las = 2)
  mtext("LOD", side = 4, line = 2, at = lodFloor + 4)
  for(diet in c("CD", "HF")) {
    Bmean <- c()
    Bstd <- c()
    Dmean <- c()
    Dstd <- c()
    for(timepoint in seq(20, 780, 30)) {
      isDIET <- ldata[which(ldata[, "Diet"] == diet & grepl("BXD", ldata[, "StrainName"])),]
      isDIET <- isDIET[which(isDIET[, "StrainName"] %in% colnames(geno)),]
      genoS <- geno[, isDIET[, "StrainName"]]

      Y <- isDIET[, "AgeAtDeath..days."]
      idx <- which(Y >= timepoint)
      Y <- Y[idx]

      strains <- isDIET[idx, "StrainName"]

      AgeAtSetUp <- as.numeric(isDIET[, "AgeAtSetUp.in.colony..days."])
      AgeAtSetUp <- as.factor(cut(AgeAtSetUp, seq(0, 1000, 180))) ### Threshold for less extreme
      AgeAtSetUp <- AgeAtSetUp[idx]

      vivarium <- as.factor(isDIET[, "vivarium"])
      vivarium <- vivarium[idx]

      null <- lmer(Y ~ AgeAtSetUp + (1 | strains), REML = FALSE)
      if(length(levels(vivarium)) >1) null <- lmer(Y ~ vivarium + AgeAtSetUp + (1 | strains), REML = FALSE)
      adj <- residuals(null, na.action=na.exclude) + mean(Y)

      gts <- as.numeric(factor(as.character(geno[ii, strains]), levels = c("B", "H", "D"))) - 2
      cat(ii, " ", timepoint, " ", diet, " ", length(adj), "==", length(gts), "\n")
      Bmean <- c(Bmean, -(mean(Y) - mean(adj[which(gts == -1)])))
      Bstd <- c(Bstd, std(adj[which(gts == -1)]))
      Dmean <- c(Dmean, -(mean(Y) - mean(adj[which(gts == 1)])))
      Dstd <- c(Dstd, std(adj[which(gts == 1)]))
    }

    if(diet != "NAM") polygon(c(seq(20, 780, 30), rev(seq(20, 780, 30))), c(Bmean + Bstd, rev(Bmean - Bstd)), col = rgb(0,0,0,0.1), border = NA)
    if(diet != "NAM") polygon(c(seq(20, 780, 30), rev(seq(20, 780, 30))), c(Dmean + Dstd, rev(Dmean - Dstd)), col = rgb(210/255,105/255,30/255,0.1), border = NA)

    points(seq(20, 780, 30), Bmean, t = "p", pch = 15+which(c("HF", "CD", "NAM") == diet), cex=1.5, col = "white")
    points(seq(20, 780, 30), Dmean, t = "p", pch = 15+which(c("HF", "CD", "NAM") == diet), col = "white", cex = 1.5)

    points(seq(20, 780, 30), Bmean, t = "b", pch = 15+which(c("HF", "CD", "NAM") == diet))
    points(seq(20, 780, 30), Dmean, t = "b", col = "chocolate", pch = 15+which(c("HF", "CD", "NAM") == diet))

    lod <- as.numeric(lods[[diet]][ii, as.character(seq(20, 780, 30))])
    points(seq(20, 780, 30), lod + lodFloor, t = "l", lwd = 2, lty = 1 + (diet == "CD"))
  }
  legend("topleft", c("B", "D"), col = c("black", "chocolate"), lwd = 2, bty = "n")
  legend("topright", c("CD", "HF"), pch = c(17, 16), bty = "n")
  legend(20, -20, c("CD", "HF"), lty = c(2, 1), lwd = 2, bty = "n", horiz = TRUE, xjust = 0, yjust = 0)

  dev.off()
}





