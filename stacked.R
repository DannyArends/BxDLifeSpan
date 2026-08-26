#
# stacked.R
#


library(RColorBrewer)
source("shared.R")

res.main <- read.table("output/HF_wMeans_Progressive.txt", sep = "\t")
res.main <- t(res.main[, -c(1:4)])

res.int <- read.table("output/CD_wMeans_Progressive.txt", sep = "\t")
res.int <- t(res.int[, -c(1:4)])

subset <- map[which(map[,1] %in% c(1:19, "X")),]
subset <- cbind(subset, cpos = NA)
gap <- 30000000
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in c(1:19, "X")){
  cl <- 1e6 * max(as.numeric(subset[which(subset[,1] == chr),4]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- (1e6 * as.numeric(subset[which(subset[,1] == chr), 4])) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

pdf("output/BeniStackedPlot.pdf", width = 14, height = 24)

op <- par(mfrow = c(2,1))
off <- 3

plot(x = c(0, cp/14.9), y = c(0, 3 + (8 * off)), t = "n", xaxt = "n", yaxt = "n", xlab = "Chromosome", ylab = "")
yo <- 8*off
abline(h = yo + off, lty = 2)
for(x in rev((1:nrow(res.main))[c(TRUE, FALSE, FALSE)])){
  i <- 1
  for(chr in c(1:19, "X")){
    onChr <- subset[which(subset[,1] == chr), 2]
    xp <- subset[onChr, "cpos"]
    xp <- c(xp[1], xp, xp[length(xp)])
    yp <- c(0, as.numeric(res.main[x, onChr]), 0)

    polygon(x = xp, y = yo + yp, 
            col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(100, 100, 0, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
            border = "white", lwd=2)
    i <- i + 1
  }
  abline(h = yo, lty = 2)
  yo <- yo - off
}



axis(1, at = chr.mids,  c(1:19, "X"))
axis(2, at = seq(.5 * off, (9 * off), off),  gsub("X", "", rownames(res.main)[(1:nrow(res.main))[c(TRUE, FALSE, FALSE)]]), las=2)
axis(3, at = chr.mids,  c(1:19, "X"))


plot(x = c(0, cp/14.9), y = c(0, 3 + (8 * off)), t = "n", xaxt = "n", yaxt = "n", xlab = "Chromosome", ylab = "")
yo <- 8*off
abline(h = yo + off, lty = 2)
for(x in rev((1:nrow(res.int))[c(TRUE, FALSE, FALSE)])){
  i <- 1
  for(chr in c(1:19, "X")){
    onChr <- subset[which(subset[,1] == chr), 2]
    xp <- subset[onChr, "cpos"]
    xp <- c(xp[1], xp, xp[length(xp)])
    yp <- c(0, as.numeric(res.int[x, onChr]), 0)

    polygon(x = xp, y = yo + yp, 
            col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(100, 100, 0, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
            border = "white", lwd=2)
    i <- i + 1
  }
  abline(h = yo, lty = 2)
  yo <- yo - off
}



axis(1, at = chr.mids,  c(1:19, "X"))
axis(2, at = seq(.5 * off, (9 * off), off),  gsub("X", "", rownames(res.main)[(1:nrow(res.main))[c(TRUE, FALSE, FALSE)]]), las=2)
axis(3, at = chr.mids,  c(1:19, "X"))
dev.off()



