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
BWs <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names=FALSE)

mm <- c()   # correlation
mls <- c()  # mean lifespan
mn <- c()   # sample size
for(diet in c("CD", "HF")) {
  mr <- c()
  mnr <- c()
  mlsr <- c()
  for(tp in colnames(BWs)[1:751]){
    iix <- which(ldata[, "AgeAtSetUp.in.colony..days."] <= as.numeric(tp) & 
                 ldata[, "AgeAtDeath..days."] >= as.numeric(tp) & 
                 grepl("BXD", ldata[, "StrainName"]) &
                 ldata[, "Diet"] == diet)
    lspan <- ldata[iix,]

    lspan <- lspan[which(rownames(lspan) %in% rownames(BWs)),]
    bspan <- BWs[which(rownames(BWs) %in% rownames(lspan)),]

    cc <- cor(bspan[rownames(lspan), as.character(tp)], lspan[rownames(lspan),"AgeAtDeath..days."], use = "pair")
    cat("", tp, ",", diet, ",", cc, "\n")
    mr <- c(mr, cc)
    mlsr <- c(mlsr, mean(lspan[rownames(lspan),"AgeAtDeath..days."]))
    mnr <- c(mnr, length(rownames(lspan)))
  }
  mm <- rbind(mm, mr)
  mn <- rbind(mn, mnr)
  mls <- rbind(mls, mlsr)
}
colnames(mm) <- colnames(BWs)[1:751]
colnames(mn) <- colnames(BWs)[1:751]
colnames(mls) <- colnames(BWs)[1:751]

mm <- mm[, which(apply(mn,2,min) > 30)]
rownames(mm) <- c("CD", "HF")

# Write as PDF for panel 1
#pdf("")
op <- par(mar = c(5,4,2,4))
plot(c(0,800), c(-0.5, 0.5), t = "n", main = "Correlation",xlab = "Days", ylab = "Correlation")
abline(h = -0.05, col = "gray")
abline(h = 0.05, col = "gray")
abline(h = 0.0)

points(as.numeric(colnames(mm)), mm[1,], t = "l")
points(as.numeric(colnames(mm)), mm[2,], t = "l", col = "blue")

points(as.numeric(colnames(mls)), (mls[1,]-mls[2,])/ 250, t = "l", lty = 2)
axis(4, at = seq(0, 0.4, 0.1), seq(0, 0.4, 0.1) * 250, las=2)
legend("topleft", c("CD", "HF"), col = c("black", "blue"), lty = 1)

