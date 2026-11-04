source("shared.R")
BWs <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)

### NAM lifespan: same filters as lifespan.R lines 13-16, 24-29 (strain means) and 55-57 (individual)
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
ldata <- ldata[which(!is.na(cut(ldata[, "AgeAtSetUp.in.colony..days."], seq(0, 1000, 180)))),]
isNAM <- ldata[which(ldata[, "Diet"] == "NAM"),]
mtable <- table(isNAM[, "StrainName"])
write.table(rownames(isNAM)[which(isNAM[, "StrainName"] %in% grep("BXD", names(mtable)[which(mtable >= 3)], value = TRUE))], "output/EarTags_NAM_lifespan_wMeans.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(isNAM)[which(grepl("BXD", isNAM[, "StrainName"]) & isNAM[, "StrainName"] %in% colnames(geno))], "output/EarTags_NAM_lifespan_ind.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

### Bodyweight: same filters as bodyweight.R lines 92-106, animals used at any mapped timepoint
bdata <- bdata[grep("^BXD", bdata[, "StrainName"]),][rownames(BWs),]
tpsA <- intersect(names(which(colSums(!is.na(BWs)) > 400)), as.character(seq(20, 780, 30)))
for(diet in c("CD", "HF", "NAM")){
  tags <- c()
  for(tp in tpsA){
    isDIET <- bdata[which(bdata[, "Diet"] == diet & !is.na(BWs[, tp])),]
    mtable <- table(isDIET[, "StrainName"])
    tags <- union(tags, rownames(isDIET)[which(isDIET[, "StrainName"] %in% names(mtable)[which(mtable >= 3)])])
  }
  write.table(tags, paste0("output/EarTags_BW_", diet, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
