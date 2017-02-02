library(tidyverse)
library(jagsUI)

rm(list = ls())

source('/home/ben/linkedModels/movementModel/moveModelWithData/moveModFunctions.R')

moveDir <- "/home/ben/linkedModels/movementModel/moveModelWithData/WB/"
setwd(moveDir)

# get original data
load(file='/home/ben/linkedModels/dataForModels/cdForMovementModelWB.RData')

# get movement model output
load( file=paste0(moveDir,"runInfo.RData") )
speciesIn <- runInfo$speciesIn
load( file=paste0(moveDir,"moveModOut_", speciesIn, ".RData") )

nItersToRun <- 5

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
cdWBPred <- cdWB %>% filter( species %in% runInfo$speciesIn, cohort %in% runInfo$cohortsIn ) %>%
                     select( tagIndex,sampleIndex,season,species,riverOrdered,enc,fOcc,lOcc,maxSampNumOB )# %>%
                  #   .[1:1000,]


iters <- sort(sample((1:out$mcmc.info$n.samples),nItersToRun))
cdWBPredByIter <- getRiverPred( cdWBPred,iters )

# check transitions
cdWBPredByIter <- cdWBPredByIter %>% group_by(tagIndex) %>% mutate(lagR = lead(riverPred),trans=paste0(riverPred,lagR))
table(cdWBPredByIter$trans)

save(cdWBPredByIter,file=paste0(moveDir,"moveModOutwPred_", speciesIn, ".RData"))

tmp <- 
cdWBPredByIter %>%
  group_by(iter) %>%
  nest()

#head(data.frame(cdWBPred[110:11631,c('tag',"tagIndex","sampleIndex",'season',"river","riverOrdered","riverPred","enc","lOcc","encObsIsOB","maxSampNumOB","pTrans","trans","flag")]),20)
#which(cdWBPred$trans==14)

