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
  write.table(res, file = paste0("output/BW_", diet, "_wMeans_Progressive.txt"), sep = "\t", quote = FALSE)
}


### Plot (copy/pasted from the plots.R)
res.main <- read.table("output/BW_CD_wMeans_Progressive.txt", sep = "\t")
res.main <- res.main[, -c(1:4)]
res.int <- read.table("output/BW_HF_wMeans_Progressive.txt", sep = "\t")
res.int <- res.int[, -c(1:4)]

library(RColorBrewer)
threshold <- 2
colz.c <- c("white", "gray", brewer.pal(9, "Greens")[-c(1:6)])

map[, "Mb"] <- 1e6 * map[, "Mb"]

chrs <- c(1:19, "X")
gap <- 10000000
chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Mb"]) }))
names(chrs.length) <- chrs
l.x <- sum(chrs.length) + gap * (length(chrs)-1)

threshold <- 1

colz.c <- c("white", "#abd9e9", "#fdb863", "#e66101", "#d01c8b", "#5e3c99")

thresholds <- c(0, 1, 2, 3, 4, 5, 100)

pdf(paste0("output/actuaryQTL_BW.pdf"), width = 24, height = 18)
  lods.c.All <- t(res.main)
  op <- par(mfrow=c(2,1))
  plot(c(1, l.x), c(1, 1+ncol(res.main)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", main = "BW CD")
  abline(h = 1:nrow(lods.c.All), col = rgb(0.5,0.5,0.5))
  chr.s <- 0
  chr.pe <- 0
  c.i <- 1
  for(chr in chrs){
    rect(chr.s, 0, chr.s + chrs.length[chr], 1+nrow(lods.c.All), col = c(rgb(1,1,1,0.5), rgb(1,1,1,0.5))[1 + (c.i %% 2)], border=NA);
    c.i <- c.i + 1;
    axis(1, at = chr.s + chrs.length[chr] / 2, chr)
    rect(chr.pe, 0, chr.s, 1+nrow(lods.c.All), col = "white", border=NA)
    abline(v = chr.s)
    abline(v = chr.s + chrs.length[chr])
    rr <- rownames(map)[which(map[, "Chr"] == chr)]
    rr.pos <- map[rr, "Mb"] + chr.s
    
    for(tp in 1:nrow(lods.c.All)){
      ii <- which(lods.c.All[tp,rr] > threshold)
      for(i in ii){
        if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
        if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
        lod <- as.numeric(lods.c.All[tp,rr][i])
        #cat(x.s," ",  x.e, " ", lod, "\n")
        cN <- max(which(thresholds < lod))
        rect(x.s, tp+0, x.e, tp+1, col = colz.c[cN], border = NA)
      }
    }
    chr.pe <- chr.s + chrs.length[chr]
    chr.s <- chr.s + chrs.length[chr] + gap
  }
  axis(2, at = 0.5 + 1:nrow(lods.c.All), gsub("X", "", rownames(lods.c.All)), las=2)
  legend("topright", fill = c(colz.c), c("0 - 2", "2 - 3", "3 - 4", "4 - 5", ">5"), bg = "white", title = "-log10P")

  lods.c.All <- t(res.int)
  plot(c(1, l.x), c(1, 1+ncol(res.main)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", main = "BW HF")
  abline(h = 1:nrow(lods.c.All), col = rgb(0.5,0.5,0.5))
  chr.s <- 0
  chr.pe <- 0
  c.i <- 1
  for(chr in chrs){
    rect(chr.s, 0, chr.s + chrs.length[chr], 1+nrow(lods.c.All), col = c(rgb(1,1,1,0.5), rgb(1,1,1,0.5))[1 + (c.i %% 2)], border=NA);
    c.i <- c.i + 1;
    axis(1, at = chr.s + chrs.length[chr] / 2, chr)
    rect(chr.pe, 0, chr.s, 1+nrow(lods.c.All), col = "white", border=NA)
    abline(v = chr.s)
    abline(v = chr.s + chrs.length[chr])
    rr <- rownames(map)[which(map[, "Chr"] == chr)]
    rr.pos <- map[rr, "Mb"] + chr.s
    
    for(tp in 1:nrow(lods.c.All)){
      ii <- which(lods.c.All[tp,rr] > threshold)
      for(i in ii){
        if(i == 1){ x.s <- chr.s; }else{ x.s <- (rr.pos[i] + rr.pos[i-1])/2; }
        if(i == length(rr)){ x.e <- chr.s + chrs.length[chr]; }else{ x.e <- (rr.pos[i] + rr.pos[i+1])/2; }
        lod <- as.numeric(lods.c.All[tp,rr][i])
        #cat(x.s," ",  x.e, " ", lod, "\n")
        cN <- max(which(thresholds < lod))
        rect(x.s, tp+0, x.e, tp+1, col = colz.c[cN], border = NA)
      }
    }
    chr.pe <- chr.s + chrs.length[chr]
    chr.s <- chr.s + chrs.length[chr] + gap
  }
  axis(2, at = 0.5 + 1:nrow(lods.c.All), gsub("X", "", rownames(lods.c.All)), las=2)
  legend("topright", fill = c(colz.c), c("0 - 2", "2 - 3", "3 - 4", "4 - 5", ">5"), bg = "white", title = "-log10P")

dev.off()



