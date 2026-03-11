#
# NAM.R - Visualize the NAM mapping
#

library(RColorBrewer)
source("shared.R")

res.nam <- read.table("output/NAM_wMeans.txt", sep = "\t")


subset <- res.nam[which(res.nam[,1] %in% c(1:19, "X")),]
subset[,"Chr"] <- as.character(subset[,"Chr"])
subset <- cbind(subset, cpos = NA)
gap <- 30
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in as.character(c(1:19, "X"))){
  cl <- max(as.numeric(subset[which(subset[,"Chr"] == chr), "Mb"]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), "Mb"]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

pdf(paste0("output/lifespan_NAM.pdf"), width = 24, height = 6)
plot(x = c(0, 1e6 * cp), y = c(0, max(res.nam [,5])), t = "n", xaxt = "n", xlab = "Chromosome", ylab = "", main = "NAM Lifespan ")
i <- 1
for(chr in as.character(c(1:19, "X"))){
  onChr <- rownames(subset[which(subset[,"Chr"] == chr),])
  xp <- 1e6 * subset[onChr, "cpos"]
  xp <- c(xp[1], xp, xp[length(xp)])
  yp <- c(0, subset[onChr,5], 0)

  polygon(x = xp, y = yp, 
          col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(0, 0, 100, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
          border = "black", lwd=1)
  i <- i + 1
}
axis(1, at = 1e6*chr.mids,  c(1:19, "X"))
dev.off()


##### interactions

res.int <- read.table("output/NAM_Interaction.txt", sep = "\t")

subset <- res.int[which(res.int[,1] %in% c(1:19, "X")),]
subset[,"Chr"] <- as.character(subset[,"Chr"])
subset <- cbind(subset, cpos = NA)
gap <- 30
chr.start <- c(0)
chr.mids <- c()
cp <- 0
for(chr in as.character(c(1:19, "X"))){
  cl <- max(as.numeric(subset[which(subset[,"Chr"] == chr), "Mb"]))
  chr.start <- c(chr.start, cl + cp + gap)
  subset[which(subset[,1] == chr), "cpos"] <- as.numeric(subset[which(subset[,1] == chr), "Mb"]) + cp
  if(chr == "X"){ chr <- 20 }
  cat(chr, " ", cl, " ", cp, " ", gap, " ", chr.start[chr], "\n")
  chr.mids <- c(chr.mids, chr.start[as.numeric(chr)] + cl/2)
  cp = cl + cp + gap
}

pdf(paste0("output/interaction_NAM.pdf"), width = 24, height = 6)
plot(x = c(0, 1e6 * cp), y = c(0, max(res.int [,5])), t = "n", xaxt = "n", xlab = "Chromosome", ylab = "", main = "NAM interaction")
i <- 1
for(chr in as.character(c(1:19, "X"))){
  onChr <- rownames(subset[which(subset[,"Chr"] == chr),])
  xp <- 1e6 * subset[onChr, "cpos"]
  xp <- c(xp[1], xp, xp[length(xp)])
  yp <- c(0, subset[onChr,5], 0)

  polygon(x = xp, y = yp, 
          col = c(rgb(0, 100, 0, 100, maxColorValue = 255), rgb(0, 0, 100, 100, maxColorValue = 255))[(i %% 2 == 0) + 1], 
          border = "black", lwd=1)
  i <- i + 1
}
axis(1, at = 1e6*chr.mids,  c(1:19, "X"))
dev.off()
