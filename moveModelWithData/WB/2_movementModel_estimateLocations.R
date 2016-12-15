library(tidyverse)
library(jagsUI)

rm(list = ls())

moveDir <- "/home/ben/movementModel/moveModelWithData/WB/"
setwd(moveDir)

# get original data
load(file='/home/ben/dataForModels/cdForMovementModelWB.RData')

# get movement model output
speciesIn <- 'bkt'
load( file=paste0(moveDir,"moveModOut_", speciesIn, ".RData") )

# Fill in missing observations of location for all rows in cdWB
#
#

cdWB <- cdWB %>% 
  group_by(tagIndex) %>%
  mutate( minOcc = min(sampleIndex), fOcc = (sampleIndex==minOcc)*1,
          maxOcc = max(sampleIndex), lOcc = (sampleIndex==maxOcc)*1
  )

cdWB <- cdWB %>% filter( species %in% runInfo$speciesIn, cohort %in% runInfo$cohortsIn )

iter <- 1

trans <- array(0,c(runInfo$nRivers,runInfo$nRivers))

cdWB$riverPred <- NA
cdWB$riverPred[1] <- as.numeric( cdWB$riverOrdered[1] )

for ( i in 1:(nrow(cdWB)-1) ){
  
  # if encountered next occ, just get river from next occ 
  if( cdWB$enc[i + 1] == 1 ){
    cdWB$riverPred[i + 1] <- as.numeric( cdWB$riverOrdered[i + 1] )
  }
  # if not encountered next occ and not the last occ this occ, then can transition
  else if( cdWB$enc[i + 1] == 0 & cdWB$lOcc[i] == 0 ){
    cdWB$riverPred[i + 1] <- which(rmultinom( 1, 1, prob = transProb[ iter,cdWB$riverPred[i],,cdWB$season[i] ] ) == 1)
    cdWB$pTrans[i] <- transProb[ iter,cdWB$riverPred[i],cdWB$riverPred[i+1],cdWB$season[i] ] # for checking
  }

}

cdWB <- cdWB %>% group_by(tag) %>% mutate(lagR = lead(riverPred),trans=paste0(riverPred,lagR))
table(cdWB$trans)
head(data.frame(cdWB[430:580,c('tag','season',"river","riverOrdered","riverPred","enc","lOcc","pTrans","trans")]),100)

