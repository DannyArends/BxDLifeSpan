### Adjusted lifespan per Diet for vivarium + AgeAtSetUp
ldata <- cbind(ldata, "Yadj" = NA)
ii <- which(ldata[, "Diet"] == "CD")
vivarium <- ldata[ii,"vivarium"]
AgeAtSetUp <- ldata[ii, "AgeAtSetUp.in.colony..days."]
AgeAtSetUp <- as.factor(cut(as.numeric(AgeAtSetUp), seq(0, 1000, 180))) ### Threshold for less extreme
Y1 <- ldata[ii,"AgeAtDeath..days."]
m1 <- lm(Y1 ~ vivarium + AgeAtSetUp, na.action = na.exclude)
i <- 1
for(x in ii){
  ldata[x, "Yadj"] <- round(as.numeric(residuals(m1) + mean(Y1))[i])
  i <- i+1
}
ii <- which(ldata[, "Diet"] == "HF")
vivarium <- ldata[ii,"vivarium"]
AgeAtSetUp <- ldata[ii, "AgeAtSetUp.in.colony..days."]
AgeAtSetUp <- as.factor(cut(as.numeric(AgeAtSetUp), seq(0, 1000, 180))) ### Threshold for less extreme
Y2 <- ldata[ii,"AgeAtDeath..days."]
m2 <- lm(Y2 ~ vivarium + AgeAtSetUp, na.action = na.exclude)
i <- 1
for(x in ii){
  ldata[x, "Yadj"] <- round(as.numeric(residuals(m2) + mean(Y2))[i])
  i <- i+1
}

BWsAdj <- matrix(NA, nrow(BWs), ncol(BWs), dimnames=list(rownames(BWs), colnames(BWs)))
ourTPS <- seq(20, 780, 30)
TPStodo <- colnames(BWs)[24:751][which(colnames(BWs)[24:751] %in% as.character(ourTPS))]

for(z in TPStodo){
  ii <- which(ldata[rownames(BWs),"Diet"] == "CD")
  vivarium <- ldata[ii,"vivarium"]
  AgeAtSetUp <- ldata[ii, "AgeAtSetUp.in.colony..days."]
  AgeAtSetUp <- as.factor(cut(as.numeric(AgeAtSetUp), seq(0, 1000, 180))) ### Threshold for less extreme
  Y <- BWs[ii, z]
  tryCatch(
  {
    m <- lm(Y ~ vivarium + AgeAtSetUp, na.action = na.exclude)
    i <- 1
    for(x in ii){
      BWsAdj[x, z] <- round(as.numeric(residuals(m) + mean(Y,na.rm = TRUE))[i])
      i <- i+1
    }
  }, error = function(msg){ })
  cat("Done",z,"\n")
}

for(z in TPStodo){
  ii <- which(ldata[rownames(BWs),"Diet"] == "HF")
  vivarium <- ldata[ii,"vivarium"]
  AgeAtSetUp <- ldata[ii, "AgeAtSetUp.in.colony..days."]
  AgeAtSetUp <- as.factor(cut(as.numeric(AgeAtSetUp), seq(0, 1000, 180))) ### Threshold for less extreme
  Y <- BWs[ii, z]
  tryCatch(
  {
    m <- lm(Y ~ vivarium + AgeAtSetUp, na.action = na.exclude)
    i <- 1
    for(x in ii){
      BWsAdj[x, z] <- round(as.numeric(residuals(m) + mean(Y,na.rm = TRUE))[i])
      i <- i+1
    }
  }, error = function(msg){ })
  cat("Done",z,"\n")
}


