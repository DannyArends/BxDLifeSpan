#
# lifespan.R - Actuarial mapping Strain means as well as individual level data mapping
# 3 Month steps (~ 100 days) hit the bodyweight timepoints
# AgeAtSetUp.in.colony..days. = Day at which the diet switch occured (HFD, as well as NAM), CD didn't switch they ate chow
#
# Vivarium - https://genenetwork.org/show_trait?trait_id=10013&dataset=BXD-LongevityPublish

library(RColorBrewer)
source("shared.R")

res.main <- read.table("output/main_ind_Progressive.txt", sep = "\t")
res.main <- res.main[, -c(1:4)]
res.int <- read.table("output/int_ind_Progressive.txt", sep = "\t")
res.int <- res.int[, -c(1:4)]


map[, "Mb"] <- 1e6 * map[, "Mb"]

chrs <- c(1:19, "X")
gap <- 10000000
chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Mb"]) }))
names(chrs.length) <- chrs
l.x <- sum(chrs.length) + gap * (length(chrs)-1)

threshold <- 2

colz.c <- c("white", "gray", brewer.pal(9, "Greens")[-c(1:6)])

thresholds <- c(0, 2, 3, 4, 5, 100)

pdf(paste0("actuaryQTL_lifespan.pdf"), width = 24, height = 18)
  lods.c.All <- t(res.main)
  op <- par(mfrow=c(2,1))
  plot(c(1, l.x), c(1, 1+ncol(res.main)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", main = "G")
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
  plot(c(1, l.x), c(1, 1+ncol(res.main)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", main = "G x Diet")
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


subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- cbind(subset, cpos = NA)
gap <- 30000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in c(1:19, "X")){
  cl <- max(as.numeric(subset[which(subset[,"Chr"] == chr), "Mb"]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), "Mb"]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

off <- 3.65
lods.cM <- t(res.main)
lods.cI <- t(res.int)

for(x in (1:nrow(lods.cM))){
  pdf(paste0("output/lifespan_D",gsub("X", "", rownames(lods.cM)[x]),".pdf"), width = 24, height = 12)
  plot(x = c(0, cp), y = c(-5, 5), t = "n", xaxt = "n", xlab = "Chromosome", ylab = "", 
       main = gsub("X", "T > ", rownames(lods.cM)[x]))
  yo <- 0 #24*off
  abline(h = yo + 2, lty = 2, col = "black")
  abline(h = yo + 3, lty = 2, col = "blue")
  abline(h = -(yo + 2), lty = 2, col = "black")
  abline(h = -(yo + 3), lty = 2, col = "blue")
  i <- 1
  for(chr in c(1:19, "X")){
    onChr <- rownames(subset[which(subset[,"Chr"] == chr),])
    xp <- subset[onChr, "cpos"]
    xp <- c(xp[1], xp, xp[length(xp)])
    yp <- c(0, as.numeric(lods.cM[x, onChr]), 0)

    polygon(x = xp, y = yo + yp, 
            col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(0, 100, 0, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
            border = "black", lwd=1)
    i <- i + 1
  }

  i <- 1
  for(chr in c(1:19, "X")){
    onChr <- rownames(subset[which(subset[,"Chr"] == chr),])
    xp <- subset[onChr, "cpos"]
    xp <- c(xp[1], xp, xp[length(xp)])
    yp <- c(0, as.numeric(lods.cI[x, onChr]), 0)

    polygon(x = xp, y = -(yo + yp), 
            col = c(rgb(255, 255, 100, 100, maxColorValue = 255), rgb(255, 255, 100, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
            border = "black", lwd=1)
    i <- i + 1
  }
  abline(h = yo, lty = 2)
  axis(1, at = chr.mids,  c(1:19, "X"))
  dev.off()
}


