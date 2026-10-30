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

source("adjust.R")

plot(colnames(BWs), apply(BWs,2,function(x){sum(!is.na(x))}))
points(colnames(BWsAdj), apply(BWsAdj,2,function(x){sum(!is.na(x))}), col = "blue", pch=19)


dir.create("output/CTL", showWarnings = FALSE, recursive = TRUE)

### CTL based on strain means first (at timepoints where we have (2 x 30) ~ 60 strains available)
for(diet in c("CD", "HF")) {
  ourTPS <- seq(20, 780, 30)
  TPStodo <- colnames(BWsAdj)[24:751][which(colnames(BWsAdj)[24:751] %in% as.character(ourTPS))]
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
    lspan <- lspan[which(rownames(lspan) %in% rownames(BWsAdj)),]
    bspan <- BWsAdj[which(rownames(BWsAdj) %in% rownames(lspan)),]

  
    cc <- cor(bspan[rownames(lspan), as.character(tp)], lspan[rownames(lspan),"Yadj"], use = "pair")

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

      mat_B <- cbind(bspan[rownames(lgtspan)[B6], as.character(tp)], lgtspan[B6, "Yadj"])
      valid_B <- complete.cases(mat_B)
      strains_B <- lgtspan[B6, "StrainName"][valid_B]
      strain_counts_B <- table(strains_B)
      weights_B <- as.numeric(1 / strain_counts_B[strains_B])
      cB <- cov.wt(apply(mat_B[valid_B, ],2,rank), wt = weights_B, cor = TRUE)$cor[1,2]

      # non - weighted
      #cDold <- cor(bspan[rownames(lgtspan)[DBA], as.character(tp)], lgtspan[rownames(lgtspan)[DBA],"AgeAtDeath..days."], use = "pair")

      mat_D <- cbind(bspan[rownames(lgtspan)[DBA], as.character(tp)], lgtspan[DBA, "Yadj"])
      valid_D <- complete.cases(mat_D)
      strains_D <- lgtspan[DBA, "StrainName"][valid_D]
      strain_counts_D <- table(strains_D)
      weights_D <- as.numeric(1 / strain_counts_D[strains_D])
      cD <- cov.wt(apply(mat_D[valid_D, ],2,rank), wt = weights_D, cor = TRUE)$cor[1,2]


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
    write.table(ctlTP, file = paste0("output/CTL/", diet, "_Wcor", tp, ".txt"), sep = "\t", quote = FALSE)
    write.table(ctlTPN, file = paste0("output/CTL/", diet, "_Wn", tp, ".txt"), sep = "\t", quote = FALSE)
  }
}

for(diet in c("CD", "HF")){
  pVals <- matrix(NA, nrow(map), length(TPStodo), dimnames = list(rownames(map), TPStodo))
  for(tp in TPStodo){
    COR <- read.table(paste0("output/CTL/", diet, "_Wcor", tp, ".txt"))
    INDs <- read.table(paste0("output/CTL/", diet, "_Wn", tp, ".txt"))
    for(x in 1:nrow(COR)){
      cor <- COR[x,]
      z <- .5*log((1.0 + cor)/(1.0 - cor))
      n <- INDs[x,]
      df <- n-(length(cor)-1)
      sumOfSq <- sum(df * z^2)
      sqOfSum <- sum(df * z)
      denom <- sum(df)
      Cv <- sumOfSq - (sqOfSum^2/ denom)

      pVals[rownames(map)[as.numeric(rownames(COR)[x])], as.character(tp)] <- sign(COR[x, 2] - COR[x, 1]) * -log10(pchisq(Cv, 1, 0, FALSE))
    }
    cat(diet, " ", tp, "\n")
  }
  write.table(round(pVals,2), paste0("output/CTL_", diet, "_LODmatrix.txt"), sep = "\t", quote = FALSE)
}

for(diet in c("CD", "HF")){
  nMat <- matrix(NA, nrow(map), length(TPStodo), dimnames = list(rownames(map), TPStodo))
  for(tp in TPStodo){
    INDs <- read.table(paste0("output/CTL/", diet, "_Wn", tp, ".txt"))
    nMat[as.numeric(rownames(INDs)), tp] <- paste(INDs[, 2], INDs[, 1], sep = ",")
  }
  write.table(nMat, paste0("output/CTL_", diet, "_Nmatrix.txt"), sep = "\t", quote = FALSE)
}

