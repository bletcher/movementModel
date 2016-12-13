library(tidyverse)
library(jagsUI)

rm(list = ls())

## defining data section ##

moveDir <- "/home/ben/movementModel/simulation/"
setwd(moveDir)

nInd <- 200
nOcc <- 10
nRivers <- 3
nSeasons <- 4

#proportion of NAs
pNA <- 0.25

trans <- 
  array(
    NA, c(nRivers,nRivers, nSeasons)
  )

trans[,,1] <- matrix( 
  c( 1, 0, 0,
     0, 1, 0,
     0, 0, 1
  ), nrow=nRivers, byrow=TRUE)
trans[,,2] <- matrix( 
  c( 1, 0, 0,
     0, 1, 0,
     0, 0, 1
  ), nrow=nRivers, byrow=TRUE)
trans[,,3] <- matrix( 
  c( 0.6, 0, 0,
     0,   1, 0,
     0.4, 0, 1
  ), nrow=nRivers, byrow=TRUE)
trans[,,4] <- matrix( 
  c( 1, 0, 0,
     0, 1, 0,
     0, 0, 1
  ), nrow=nRivers, byrow=TRUE)

nAdapt = 100
nIter = 1000
nBurn = 500

## end of defining data section ##

dat <- matrix(NA,nInd,nOcc)#, dimnames=list(1:nOcc, 1:T)) 
transPairs <- matrix(NA,nInd,nOcc)

for (i in 1:nInd){
  dat[i,1] <- sample( 1:nRivers, 1 ) 
  for (t in 1:(nOcc-1)){

      seas <- t %% nSeasons
      if(seas == 0) seas <- nSeasons

      dat[i,t+1] <- which(rmultinom(1, 1, prob = trans[ ,dat[i,t],seas ]) == 1)
      transPairs[i,t] <- paste0(dat[i,t],dat[i,t+1])
    
  }
}

# check transitions
for (i in 1:(nOcc-1)){
  print(c(i,table(transPairs[,i])))
}


d <- gather(as.data.frame(t(dat))) # to long format
names(d) <- c("id","r")
d$occ <- rep(1:nOcc,nInd)
d$seas <- rep( rep((1:4), nOcc/nSeasons, length.out=nOcc), nInd )

d <- d %>% 
  group_by(id) %>%
  mutate( minOcc = min(occ), fOcc = (occ==minOcc)*1,
          maxOcc = max(occ), lOcc = (occ==maxOcc)*1
        )
evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)

d$rNA <- ifelse(runif(nrow(d),0,1) < pNA & d$fOcc != 1, NA, d$r) #always have a r for fOcc

sink("moveModSim.jags")
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

data <- list( riverDATA=d$rNA, 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$occ, season = d$seas,
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
            model.file = paste0(moveDir,"moveModSim.jags"),
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

simInfo <- list(nInd=nInd,nOcc=nOcc,nRivers=nRivers,nSeasons=nSeasons)
save(out,simInfo,file=paste0(moveDir,"simMoveOut.RData"))
