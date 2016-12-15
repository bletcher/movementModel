library(tidyverse)
library(jagsUI)

rm(list = ls())
source('/home/ben/movementModel/moveModelWithData/moveModFunctions.R')

## defining data section ##

moveDir <- "/home/ben/movementModel/moveModelWithData/WB/"
setwd(moveDir)

load(file='/home/ben/dataForModels/cdForMovementModelWB.RData')

###
# Variables to change
speciesIn <- 'bkt'
cohortsIn <- c(2002) #sort(unique(cdWB$cohort)) # for all cohorts

nAdapt = 100
nIter = 1000
nBurn = 500
nThin = 5

###
#prepare and filter the data
cdWB <- cdWB %>% filter( knownZ == 1 )

#limit data to more than one obs per fish
counts <- cdWB %>% group_by(tag) %>% summarize(n=n()) %>% arrange(desc(n))
cdWB <- left_join(cdWB,counts) %>%
        filter(n > 1)  # remove single observations 

d <- cdWB %>% 
  group_by(tagIndex) %>%
  mutate( minOcc = min(sampleIndex), fOcc = (sampleIndex==minOcc)*1,
          maxOcc = max(sampleIndex), lOcc = (sampleIndex==maxOcc)*1,
          lagR = lead(as.numeric(riverOrdered)),trans=paste0(as.numeric(riverOrdered),lagR)
        )

d <- d %>% filter( species %in% speciesIn, cohort %in% cohortsIn )

evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)

nInd <- length(unique(d$tag))
nOcc <- length(unique(d$sampleIndex))
nRivers <- length(unique(d$river))
nSeasons <- length(unique(d$season))

###
# prepare for model run
data <- list( riverDATA=as.numeric(d$riverOrdered), 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$sampleIndex, season = d$season,
              nEvalRows=nEvalRows, evalRows = evalRows,
              nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
              nSeasons = nSeasons)

psiBetaInits <- array(rnorm(nRivers^2*nSeasons,0,2.25),c(nRivers,nRivers,nSeasons))
psiBetaInits[1:3,4,] <- NA # need init value of NA for fixed parameters

inits <- function(){
  list(psiBeta=psiBetaInits)
}

params <- c("psiBeta01")

out <- jags(data = data,
            inits = inits,
            parameters.to.save = params,
            model.file = paste0(moveDir,"moveMod.jags"),
            n.chains = 3,
            n.adapt = nAdapt,
            n.iter = nIter,
            n.burnin = nBurn,
            n.thin = nThin,
            parallel = TRUE
           )
###
# str(out[1]$sims.list$riverDATA[,]) gives num [1:(niter+n.burn)/n.thin, 1:nEvalRows] for chain # [1]
# get all iterations for the 1st fish's first obs as out[1]$sims.list$riverDATA[,(1-1)*nOcc+1]
# get all iterations for the 2nd fish's first obs as out[1]$sims.list$riverDATA[,(2-1)*nOcc+1]
#
# to leave out the burnin iters for the 1st fish's first obs as 
#    out[1]$sims.list$riverDATA[((out$mcmc.info$n.burn/out$mcmc.info$n.thin)+1):((out$mcmc.info$n.burn+out$mcmc.info$n.iter)/out$mcmc.info$n.thin),(1-1)*nOcc+1]
#
# to get all evalRows for 1 iter for non burn in iters:
#    out[1]$sims.list$riverDATA[((out$mcmc.info$n.burn/out$mcmc.info$n.thin)+1),]

###
# index for psiBeta01: [1:nIters, 1:nrivers, 1:nRivers, 1:nSeasons]
# to get all iters for psiBeta01: out[1]$sims.list$psiBeta01
# nIters = (nIter+nBurnin)/nThin
#
# there is a non-zero prob of transitioning based on priors (see above), even for impossible transistions (WB->OB), 
# so provide a cuttoff (0.005, see hist below) below which probs are 0
# hist(tmp[,4,1,])
rm(transProb)
transProb <- fillPsiBeta01(1,0.0075,data,out)

runInfo <- list(nInd=nInd,nOcc=nOcc,nRivers=nRivers,nSeasons=nSeasons, speciesIn=speciesIn, cohortsIn=cohortsIn )
save(d,data,out,transProb,runInfo,file=paste0(moveDir,"moveModOut_", speciesIn, ".RData"))

