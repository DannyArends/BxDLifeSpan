#
# Shared.R - Code used to load and preprocess our data

#setwd("C:/Github/BxDLifespan/")
ldata <- read.table("data/BxD_Lifespan.txt", sep = "\t", header = TRUE, row.names=2)
bdata <- read.table("data/BxD_Bodyweights.txt", sep = "\t", header = TRUE)

### One eartag is found attached to two different mice
eartag <- names(which(table(bdata[,1]) > 1))
bdata <- bdata[-which(bdata[,1] == eartag),]
ldata <- ldata[-which(rownames(ldata) == eartag),]

### Only take animals which are not "Harvested"
iix <- which(bdata[, "Status"] %in% c("Harvest", "Other"))
bdata <- bdata[-iix,]

### Still some missing mice (Yes bw, no lifespan)
rownames(bdata) <- bdata[,1]

### Mice which are matched
bdata <- bdata[which(rownames(bdata) %in% rownames(ldata)),]
ldata <- ldata[which(rownames(ldata) %in% rownames(bdata)),]

### Sort based on the order in lifespan
bdata <- bdata[rownames(ldata),]

### Write out the animals we used:
cat(bdata[, "EarTagNumberCurrent"], sep = "\n", file = "output/EarTagNumberCurrent.txt")


### Load in genotypes
geno <- read.table("http://files.genenetwork.org/current/GN600/BXD_current_rev050423.geno", 
                   sep = "\t", skip=23, header = TRUE, na.strings = c("U", "", "NA"))
map <- geno[, 1:4]
geno <- geno[, -c(1:4)]
rownames(geno) <- map[,2]
rownames(map) <- map[,2]

