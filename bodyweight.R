#
# bodyweight.R Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
# AgeAtSetUp == DateWeightBaseline == DateDietStart 
# Weight00.Baseline ~~ AgeAtSetUp


library(lme4)
library(parallel)

source("shared.R")

bdata <- bdata[grep("^BXD", bdata[, "StrainName"]),]

### Growth curves and figuring out timepoints and good samples to keep
plot(c(0, 1200), c(0, 120), t = "n")
alltps <- c()
good <- c()
for(x in sample(rownames(bdata))) {
  tps <- c("AgeAtSetUp", colnames(bdata)[grep("AgeAtWeight", colnames(bdata))])
  tp <- as.numeric(bdata[x,tps])
  wps <- colnames(bdata)[grep("^Weight", colnames(bdata))]
  weight <- as.numeric(bdata[x,wps])

  weights.na <- which(is.na(weight))

  alltps <- unique(c(alltps, tp[-weights.na]))

  if(length(tp[-weights.na]) > 2) good <- c(good, x)

  points(tp[-weights.na], weight[-weights.na], t = "l", col = c(rgb(1,0,1,0.5), rgb(0,1,1,0.5))[c("CD", "HF") == bdata[x, "Diet"] ])
}


### Create a unified bodyweight matrix
alltps <- sort(alltps)
BWs <- matrix(NA, length(good), length(alltps), dimnames = list(as.character(good), as.character(alltps)))
for(x in good){
  tps <- c("AgeAtSetUp", colnames(bdata)[grep("AgeAtWeight", colnames(bdata))])
  tp <- as.numeric(bdata[x,tps])
  wps <- colnames(bdata)[grep("^Weight", colnames(bdata))]
  weight <- as.numeric(bdata[x,wps])

  weights.na <- which(is.na(weight))
  BWs[x, as.character(tp[-weights.na])] <-  weight[-weights.na]
}


### Linear interpolation on BWs
for(x in good){
  tps <- colnames(BWs)[which(!is.na(BWs[x,]))]
  y <- BWs[x, tps]
  tps <- as.numeric(tps)
  tp.min <- min(tps)
  tp.max <- max(tps)

  xs <- which(colnames(BWs) == tp.min)
  xe <- which(colnames(BWs) == tp.max)

  res.a <- approx(tps, y, as.numeric(colnames(BWs)[xs:xe]))
  plot(c(0, 1200), c(0, 120), t = "n")
  points(res.a$x, res.a$y); points(tps, y, col = "red", pch = 16)
  BWs[x, as.character(res.a$x)] <- round(res.a$y,1)
  cat(x, "[",length(which(!is.na(BWs[x,]))),"]: ", tp.min,"-", tp.max, "\n")
}

# Subset to only the BWs animals
bdata <- bdata[rownames(BWs),]

# Sample sizes
diet <- bdata[, "Diet"]
CD <- which(diet == "CD")
HF <- which(diet == "HF")
nSamples <- nrow(BWs) - apply(apply(BWs, 2, is.na),2,sum)
n.CD <- length(CD) - apply(apply(BWs[CD,], 2, is.na),2,sum)
n.HF <- length(HF) - apply(apply(BWs[HF,], 2, is.na),2,sum)

timepoints <- as.numeric(colnames(BWs))
plot(timepoints, nSamples, main = "Sample Size", pch = 19, ylab = "nSamples", xlab = "Time")
points(timepoints, n.CD, pch = 19, col = rgb(1,0,1,0.5))
points(timepoints, n.HF, pch = 19, col = rgb(0,1,1,0.5))
legend("topright", c("All", "CD", "HF"), pch = 19, col = c(1, rgb(1,0,1,0.5),  rgb(0,1,1,0.5)))

# Reasonable timepoints
mapAble <- names(which(nSamples > 400))
tpsA <- mapAble[which(mapAble %in% as.character(seq(30, 780, 30)))]

### Mapping Weighted Strain means

for(diet in c("CD", "HF")) {
  pvalsM <- c()
  op <- par(mfrow = c(2,1))
  for(tp in tpsA){
    # Animals on the diet with bodyweight available
    isDIET <- bdata[which(bdata[, "Diet"] == diet & !is.na(BWs[,tp])),]
    mtable <- table(isDIET[, "StrainName"])
    mtable <- mtable[which(mtable >= 3)]
    strains <- names(mtable)
    bxds <- strains[grep("BXD", strains)]
    weights <- mtable[bxds]

    # Compute strain means for the strains that have at least 3 animals
    strainM <- round(unlist(lapply(bxds, function(x){
      isStrain <- rownames(isDIET[which(isDIET[, "StrainName"] == x),])
      return(mean(BWs[isStrain, tp], na.rm = TRUE))
    })), 1)
    names(strainM) <- bxds

    genoS <- geno[,names(strainM)]

    # Mapping
    pvals <- c()
    for(x in c(1:nrow(geno))) {
      gts <- as.numeric(factor(as.character(genoS[x, ]), levels = c("B", "H", "D"))) - 2
      pvals <- c(pvals, anova(lm(strainM ~ gts, weights = weights))[[5]][1])
      if(x %% 1000 == 1) cat(x, "\n")
    }
    plot(-log10(pvals), col = as.numeric(as.factor(map[, "Chr"])), main = paste0(diet, "@", tp, ",n=", length(bxds)))
    pvalsM <- cbind(pvalsM, pvals)
  }
  res <- cbind(map, -log10(pvalsM))
  colnames(res) <- c("Chr","Locus","cM","Mb", tpsA)
  write.table(res, file = paste0("output/BW_", diet, "_wMeans.txt"), sep = "\t", quote = FALSE)
}

