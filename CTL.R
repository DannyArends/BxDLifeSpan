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
ldata <- cbind(ldata, AgeAtSetUpGroup = cut(ldata[, "AgeAtSetUp.in.colony..days."], seq(0, 1000, 180)))

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

### CTL based on strain means first (at timepoints where we have (2 x 30) ~ 60 strains available)
##
for(diet in c("HF")) {
  ourTPS <- seq(20, 780, 30)
  TPStodo <- colnames(BWs)[24:751][which(colnames(BWs)[24:751] %in% as.character(ourTPS))]
  mCTL <- vector("list", length(TPStodo))
  names(mCTL) <- TPStodo
  mCTLN <- vector("list", length(TPStodo))
  names(mCTLN) <- TPStodo

  for(tp in TPStodo) {
    iix <- which(ldata[, "AgeAtSetUp.in.colony..days."] <= as.numeric(tp) & 
                 ldata[, "AgeAtDeath..days."] >= as.numeric(tp) & 
                 grepl("BXD", ldata[, "StrainName"]) &
                 ldata[, "Diet"] == diet)
    lspan <- ldata[iix,]
    lspan <- lspan[which(rownames(lspan) %in% rownames(BWs)),]
    bspan <- BWs[which(rownames(BWs) %in% rownames(lspan)),]
    cc <- cor(bspan[rownames(lspan), as.character(tp)], lspan[rownames(lspan),"AgeAtDeath..days."], use = "pair")

    ctlTP <- c()
    ctlTPN <- c()
    mymarkers <- c()
    for(m in 1:nrow(geno)) {
      strains <- lspan[, "StrainName"]
      goodS <- which(strains %in% colnames(geno))
      strains <- strains[goodS]

      lgtspan <- cbind(lspan[goodS, ], gts = as.character(geno[m, strains]))
      B6 <- which(lgtspan[, "gts"] == "B")
      DBA <- which(lgtspan[, "gts"] == "D")
      if(length(B6) < 20) next;
      if(length(DBA) < 20) next;

      # non - weighted
      #cBold <- cor(bspan[rownames(lgtspan)[B6], as.character(tp)], lgtspan[rownames(lgtspan)[B6],"AgeAtDeath..days."], use = "pair")

      mat_B <- cbind(bspan[rownames(lgtspan)[B6], as.character(tp)], lgtspan[B6, "AgeAtDeath..days."])
      valid_B <- complete.cases(mat_B)
      strains_B <- lgtspan[B6, "StrainName"][valid_B]
      strain_counts_B <- table(strains_B)
      weights_B <- as.numeric(1 / strain_counts_B[strains_B])
      cB <- cov.wt(mat_B[valid_B, ], wt = weights_B, cor = TRUE)$cor[1,2]

      # non - weighted
      #cDold <- cor(bspan[rownames(lgtspan)[DBA], as.character(tp)], lgtspan[rownames(lgtspan)[DBA],"AgeAtDeath..days."], use = "pair")

      mat_D <- cbind(bspan[rownames(lgtspan)[DBA], as.character(tp)], lgtspan[DBA, "AgeAtDeath..days."])
      valid_D <- complete.cases(mat_D)
      strains_D <- lgtspan[DBA, "StrainName"][valid_D]
      strain_counts_D <- table(strains_D)
      weights_D <- as.numeric(1 / strain_counts_D[strains_D])
      cD <- cov.wt(mat_D[valid_D, ], wt = weights_D, cor = TRUE)$cor[1,2]


      #cat("", tp, ",", diet, ",", m, "=", cc, ",",cB, ",", cD, "\n")
      ctlTP <- rbind(ctlTP, c(cB, cD))
      ctlTPN <- rbind(ctlTPN, c(length(weights_B), length(weights_D)))
      mymarkers <- c(mymarkers, m)
    }
    cat("", tp, ",", diet, "DONE\n")
    rownames(ctlTP) <- mymarkers
    rownames(ctlTPN) <- mymarkers
    mCTL[[tp]] <- ctlTP
    mCTLN[[tp]] <- ctlTPN
  }
}

for(x in 1:length(mCTL)){
  write.table(mCTL[[x]], file = paste0("output/CTL/HF_Wcor", names(mCTL)[x], ".txt"), sep = "\t", quote = FALSE)
  write.table(mCTLN[[x]], file = paste0("output/CTL/HF_Wn", names(mCTLN)[x], ".txt"), sep = "\t", quote = FALSE)
}

tp <- 5

op <- par(mfrow=c(2,1))
plot(c(0, nrow(mCTL[[tp]])), c(-1, 1), t = "n", main = "170Days")
abline(h = cc)
points(mCTL[[tp]][,1], t = "l", col = "gray")
points(mCTL[[tp]][,2], t = "l", col = "red")
legend("topleft", c("B", "D"), lty=1, col = c("gray", "red"))

plot(c(0, nrow(mCTLN[[tp]])), c(0, 500), t = "n", main = "170Days")
abline(h = cc)
points(mCTLN[[tp]][,1], t = "l", col = "gray")
points(mCTLN[[tp]][,2], t = "l", col = "red")
legend("topleft", c("B", "D"), lty=1, col = c("gray", "red"))



CDn <- read.table("output/CTL/CD_cor170.txt")
CDw <- read.table("output/CTL/CD_Wcor170.txt")

op <- par(mfrow = c(1,2))
plot(c(-0.7, 0.3), c(-0.7, 0.3), t = "n", xlab = "No strain", ylab = "Weighted strain", main = "B - allele carriers")
points(CDn[,1], CDw[,1])
plot(c(-0.7, 0.3), c(-0.7, 0.3), t = "n", xlab = "No strain", ylab = "Weighted strain", main = "D - allele carriers")
points(CDn[,2], CDw[,2])



pVals <- matrix(NA, dim(map), length(TPStodo))
rownames(pVals) <- rownames(map)
colnames(pVals) <- TPStodo
for(tp in TPStodo){
  COR <- read.table(paste0("output/CTL/HF_cor",tp,".txt"))
  INDs <- read.table(paste0("output/CTL/HF_n",tp,".txt"))
  pC <- c()
  for(x in 1:nrow(COR)){
    cor <- COR[x,]
    z <- .5*log((1.0 + cor)/(1.0 - cor))
    n <- INDs[x,]
    df <- n-(length(cor)-1)
    sumOfSq <- sum(df * z^2)
    sqOfSum <- sum(df * z)
    denom <- sum(df)
    Cv <- sumOfSq - (sqOfSum^2/ denom)

    pVals[rownames(map)[as.numeric(rownames(COR)[x])], tp] <- pchisq(Cv, 1, 0, FALSE)
    #cat(x, " ", as.numeric(cor[1]), " ", as.numeric(cor[2])," ",  Cv, "\n")
  }
}

write.table(pVals, "output/CTL_HF_Pmatrix.txt", sep = "\t", quote=FALSE)

pVals <- read.table("output/CTL_HF_Pmatrix.txt", sep = "\t")

image(-log10(apply(pVals,2,as.numeric)), breaks = c(0,1,2,3,4,5,6,10), 
      col = c("white", "#E1E1E1", "#D0D0D0", "#00C0C0","#00B0B0","#00A0A0","#900000"), xaxt="n", yaxt = "n")
axis(2, at = (1:7 - 0.5) / 7, seq(50, 780, 30)[1:8], las = 2)
