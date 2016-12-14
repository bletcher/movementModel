library(tidyverse)
library(jagsUI)

rm(list = ls())

## defining data section ##

moveDir <- "/home/ben/movementModel/moveModelWithData/WB/"
setwd(moveDir)

load(file='/home/ben/dataForModels/cdForMovementModelWB.RData')

nAdapt = 100
nIter = 1000
nBurn = 500

cdWB <- cdWB %>% filter( knownZ ==1 )

counts <- cdWB %>% group_by(tag) %>% summarize(n=n()) %>% arrange(desc(n))
cdWB <- left_join(cdWB,counts) %>%
        filter(n > 1)  # remove single observations 

d <- cdWB %>% 
  group_by(tagIndex) %>%
  mutate( minOcc = min(sampleIndex), fOcc = (sampleIndex==minOcc)*1,
          maxOcc = max(sampleIndex), lOcc = (sampleIndex==maxOcc)*1
        )
d$riverIndex <- as.numeric(factor(d$river))

d <- d %>% filter(cohort == 2002 & species == 'bkt')

evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)

nInd <- length(unique(d$tag))
nOcc <- length(unique(d$sampleIndex))
nRivers <- length(unique(d$river))
nSeasons <- length(unique(d$season))

sink("moveMod.jags")
cat("
    model {
    
    # Priors and constraints

    for( t in 1:nSeasons ) {
      for(r in 1:(nRivers)){
        for(r2 in 1:(nRivers)){
          
          psiBeta[r,r2,t] ~ dnorm( 0,0.1 ) #1/2.25) 
          psiBeta01[r,r2,t] <- exp( psiBeta[r,r2,t] )/( exp( psiBeta[r,r2,t] ) + 1 )

          # Note that with precision priors on psiBeta of 1/2.25 (value from the paper) estimates of psiBeta01 when r,r2,t sets have no transitions 
          # don't bump up against 0, but they are > 0 (~0.005). Smaller priors (<0.01) give spikey estimates of psiBeta01 with poor r-hats.
          # 0.1 is a reasonable compromise. Alternavitvely could stick with 1/2.25 (which samples better) and convert estimates of < 0.005 to 0

        }
      }
    }

    ############################
    ######## psi model #########
    ############################
    # note that 'riverDATA' is equivalent to 'r1'    

    for( i in 1:nEvalRows ){

      sumPsi[evalRows[i]] <- sum(ePsi[evalRows[i],])
    
      for( r2 in 1:nRivers ){
        # normal priors on logit
        lpsi[evalRows[i],r2] <- psiBeta[ riverDATA[evalRows[i]], r2, season[evalRows[i]] ] # psiBeta[r,r,] is not estimated. Get by difference.    
    
        ePsi[evalRows[i],r2] <- exp(lpsi[evalRows[i],r2])*(1-(riverDATA[evalRows[i]]==r2))
    
        #Constrain each set of psi's to sum to one
        psi[evalRows[i],r2] <- ( ePsi[evalRows[i],r2] / ( 1 + sumPsi[evalRows[i]]) ) * ( 1 - (riverDATA[evalRows[i]]==r2) ) # prob of moving to r2
                                                + ( 1 / ( 1 + sumPsi[evalRows[i]]) ) * (      riverDATA[evalRows[i]]==r2 )  # prob of staying in current location
      }
    }

    ############################
    ##### Likelihoods ##########
    ############################

#    for(i in 1:nFirstObsRows){
#       zRiv[ firstObsRows[i] ] <- riverDATA[ firstObsRows[i] ] 
#   }
    
    for(i in 1:nEvalRows){

      # State of location (zRiv)
      riverDATA[ evalRows[i] + 1 ] ~ dcat( psi[ evalRows[i], ]  )
#      zRiv[ evalRows[i] + 1 ] <- riverDATA[ evalRows[i] + 1 ] # this was used to truncate occasions to when z==1 
    }

  } # model
    ",fill = TRUE)
sink()

data <- list( riverDATA=d$riverIndex, 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$sampleIndex, season = d$season,
              nEvalRows=nEvalRows, evalRows = evalRows,
              nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
              nSeasons = nSeasons)

inits <- function(){
  list(psiBeta=array(rnorm(nRivers^2*nSeasons,0,2.25),c(nRivers,nRivers,nSeasons)))
}

params <- c("psiBeta01", "riverDATA")

out <- jags(data = data,
            inits = inits,
            parameters.to.save = params,
            model.file = paste0(moveDir,"moveMod.jags"),
            n.chains = 3,
            n.adapt = nAdapt,
            n.iter = nIter,
            n.burnin = nBurn,
            n.thin = 5,
            parallel = TRUE
           )
# str(out[1]$sims.list$riverDATA[,]) gives num [1:(niter+n.burn)/n.thin, 1:nEvalRows] for chain # [1]
# get all iterations for the 1st fish's first obs as out[1]$sims.list$riverDATA[,(1-1)*nOcc+1]
# get all iterations for the 2nd fish's first obs as out[1]$sims.list$riverDATA[,(2-1)*nOcc+1]

# to leave out the burnin iters for the 1st fish's first obs as 
#    out[1]$sims.list$riverDATA[((out$mcmc.info$n.burn/out$mcmc.info$n.thin)+1):((out$mcmc.info$n.burn+out$mcmc.info$n.iter)/out$mcmc.info$n.thin),(1-1)*nOcc+1]

# to get all evalRows for 1 iter for non burn in iters:
#    out[1]$sims.list$riverDATA[((out$mcmc.info$n.burn/out$mcmc.info$n.thin)+1),]

runInfo <- list(nInd=nInd,nOcc=nOcc,nRivers=nRivers,nSeasons=nSeasons)
save(out,runInfo,file=paste0(moveDir,"moveModOut.RData"))
