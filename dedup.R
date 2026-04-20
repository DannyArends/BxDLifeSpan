
library(lme4)
library(parallel)

source("shared.R")

strains <- lspan[, "StrainName"]
goodS <- which(strains %in% colnames(geno))
strains <- strains[goodS]

tokeep <- c()
x <- 1

while(x < (nrow(geno)-1)){
  M <- paste0(geno[x,strains], collapse="")
  tokeep <- c(tokeep, x)
  cat(x, ", kept: ", length(tokeep), "\n")
  y <- x + 1
  while(y < nrow(geno)){
    if(paste0(geno[y,strains], collapse="") == M){
      y <- y + 1
    } else {
      tokeep <- c(tokeep, y-1)
      break;
    }
  }
  x <- y
}


geno <- geno[tokeep, ]
