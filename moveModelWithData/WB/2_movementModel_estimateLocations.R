library(tidyverse)
library(jagsUI)

rm(list = ls())

moveDir <- "/home/ben/movementModel/moveModelWithData/WB/"
setwd(moveDir)

# get original data
load(file='/home/ben/dataForModels/cdForMovementModelWB.RData')

# get movement model output
load( file=paste0(moveDir,"runInfo.RData") )
speciesIn <- runInfo$speciesIn
load( file=paste0(moveDir,"moveModOut_", speciesIn, ".RData") )

# Fill in missing observations of location for all rows in cdWB
#
#

cdWB <- cdWB %>% 
  group_by(tagIndex) %>%
  mutate( minOcc = min(sampleIndex), fOcc = (sampleIndex==minOcc)*1,
          maxOcc = max(sampleIndex), lOcc = (sampleIndex==maxOcc)*1,
        #  lastObsIsOBT = (lOcc == 1 & ( riverOrdered == 'wb obear') ),
          encObsIsOB = (enc == 1 & ( riverOrdered == 'wb obear') )
          
          
         )

# What sample# is the last obs in O'Bear for each fish?
lastOB <- cdWB %>%
  group_by(tagIndex) %>%
  summarize( #lastObsIsOB = any( lastObsIsOBT, na.rm = TRUE ),
             maxSampNumOB = max(sampleIndex[encObsIsOB])  )

cdWB <- left_join( cdWB, lastOB )

# filter to species and cohorts used in estimation
cdWB <- cdWB %>% filter( species %in% runInfo$speciesIn, cohort %in% runInfo$cohortsIn )

iter <- 1

trans <- array(0,c(runInfo$nRivers,runInfo$nRivers))

cdWB$pTrans <- NA
cdWB$riverPred <- NA
cdWB$riverPred[1] <- as.numeric( cdWB$riverOrdered[1] )

cdWB$flag <- NA

for ( i in 1:(nrow(cdWB)-1) ){
  if(i %% 1000 == 0) print(i)

  # if encountered next occ, just get river from the next occ 
  if( cdWB$enc[i + 1] == 1 ){
    cdWB$riverPred[i + 1] <- as.numeric( cdWB$riverOrdered[i + 1] )
    cdWB$flag[i] <- 1
  }
  # if not encountered next occ and not the last occ this occ, then can transition
  else if( cdWB$enc[i + 1] == 0 & cdWB$lOcc[i] == 0 ){
    cdWB$riverPred[i + 1] <- which(rmultinom( 1, 1, prob = transProb[ iter,cdWB$riverPred[i],,cdWB$season[i] ] ) == 1)
    cdWB$pTrans[i] <- transProb[ iter,cdWB$riverPred[i],cdWB$riverPred[i+1],cdWB$season[i] ] # for checking
    
    cdWB$flag[i] <- 2
    
    # fix river to OB if the last OB obs and the fish was not caught
    # maxSampNumOB is -Inf if fish never in OB, finite value is last observed OB sampleNum
    if( cdWB$enc[i + 1] == 0 & is.finite(cdWB$maxSampNumOB[i]) & cdWB$sampleIndex[i] < cdWB$maxSampNumOB[i] ) {
      cdWB$riverPred[i + 1] <- 4 
      cdWB$pTrans[i] <- 0
      
      cdWB$flag[i] <- 3
    }
  }
}

# check transitions
cdWB <- cdWB %>% group_by(tag) %>% mutate(lagR = lead(riverPred),trans=paste0(riverPred,lagR))
table(cdWB$trans)
#head(data.frame(cdWB[1:11631,c('tag',"tagIndex","sampleIndex",'season',"river","riverOrdered","enc","encObsIsOB","maxSampNumOB")]),20)

head(data.frame(cdWB[1:11631,c('tag',"tagIndex","sampleIndex",'season',"river","riverOrdered","riverPred","enc","lOcc","encObsIsOB","maxSampNumOB","pTrans","trans","flag")]),20)

which(cdWB$trans==14)

