library(RColorBrewer)
source("shared.R")

res.CD <- read.table("output/CD_wMeans_Progressive.txt", sep = "\t")
res.CD <- res.CD[, -c(1:4)]

res.HF <- read.table("output/HF_wMeans_Progressive.txt", sep = "\t")
res.HF <- res.HF[, -c(1:4)]


map[, "Mb"] <- 1e6 * map[, "Mb"]

chrs <- c(1:19, "X")
gap <- 10000000
chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x,"Mb"]) }))
names(chrs.length) <- chrs
l.x <- sum(chrs.length) + gap * (length(chrs)-1)

threshold <- 1

colz.c <- c("white", "#abd9e9", "#fdb863", "#e66101", "#d01c8b", "#5e3c99")

thresholds <- c(0, 1, 2, 3, 4, 5, 100)

pdf(paste0("output/actuaryQTL_lifespan_CD_HF.pdf"), width = 24, height = 18)
  lods.c.All <- t(res.CD)
  op <- par(mfrow=c(2,1))
  plot(c(1, l.x), c(1, 1+ncol(res.CD)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", 
       main = "CD - Weighted strain mean mapping (Main effect, n/strain >= 3)")
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
  legend("topright", fill = c(colz.c), c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", ">5"), bg = "white", title = "-log10P")

  lods.c.All <- t(res.HF)
  plot(c(1, l.x), c(1, 1+ncol(res.HF)), t = 'n', xlab = "", xaxt='n', ylab="", yaxt='n', xaxs="i", yaxs="i", 
       main = "HF - Weighted strain mean mapping (Main effect, n/strain >= 3)")
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
  legend("topright", fill = c(colz.c), c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", ">5"), bg = "white", title = "-log10P")

dev.off()

