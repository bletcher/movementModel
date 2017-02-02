
fillPsiBeta01 <- function(chainToUse,cutoff0,data,out){
  
  psiBeta01 <- out$sims.list$psiBeta01   
  
  rivs <- 1:data$nRivers
  nRivers <- data$nRivers
  nSeasons <- data$nSeasons
  nIter <- out$mcmc.info$n.iter
  nSamples <- out$mcmc.info$n.samples

  psiBeta01Filled <- array(NA,c(nSamples,nRivers,nRivers,nSeasons))
  ii=0
  for( iter in 1:nSamples ){
    for( s in 1:data$nSeasons ){
      for( r1 in 1:data$nRivers ){
        for( r2 in 1:data$nRivers ){
          ii=ii+1
          print(paste("Iter = ",iter,s,r1,r2,ii))
          if( psiBeta01[iter,r1,r2,s] < cutoff0 ) psiBeta01[iter,r1,r2,s] <- 0
          
          if (r1 == r2){
            notR <- rivs[-r1]
            notRSum <- 0
            for(r in 1:(nRivers-1)){ notRSum <- notRSum + psiBeta01[iter,r1,notR[r],s] }
            psiBeta01Filled[iter,r1,r2,s] <- 1 - notRSum
            
          }
          
          else psiBeta01Filled[iter,r1,r2,s] <- psiBeta01[iter,r1,r2,s]
          
        }
      }
    }
  }
  return(psiBeta01Filled)
}

getRiverPred <- function (cd, iters) {
  
  cdOut <- data.frame()
  
  trans <- array(0,c(runInfo$nRivers,runInfo$nRivers))
  
  cd$iter <- NA
  # cd$pTrans <- NA
  cd$riverPred <- NA
  cd$riverPred[1] <- as.numeric( cd$riverOrdered[1] )
  
  cd$flag <- NA
  
  
  for ( iter in iters ){
    startTime <- Sys.time()
    for ( i in 1:(nrow(cd)-1) ){
      if(i %% 10000 == 0) print(c(iter,i))
      
      # if encountered next occ, just get river from the next occ 
      if( cd$enc[i + 1] == 1 ){
        cd$riverPred[i + 1] <- as.numeric( cd$riverOrdered[i + 1] )
        cd$flag[i] <- 1
      }
      
      # if not encountered next occ and not the last occ this occ, then can transition
      else if( cd$lOcc[i] == 0 ){
        tP <- transProb[ iter,cd$riverPred[i],,cd$season[i] ]
        cd$riverPred[i + 1] <- which(rmultinom( 1, 1, prob = tP ) == 1)
        # cd$pTrans[i] <- transProb[ iter,cd$riverPred[i],cd$riverPred[i+1],cd$season[i] ] # for checking
        
        cd$flag[i] <- 2
        
        # fix river to OB if the last OB obs and the fish was not caught
        # maxSampNumOB is -Inf if fish never in OB, finite value is last observed OB sampleNum
        if( is.finite(cd$maxSampNumOB[i]) & cd$sampleIndex[i] < cd$maxSampNumOB[i] ) {
          cd$riverPred[i + 1] <- 4 
          #         cd$pTrans[i] <- 0
          
          cd$flag[i] <- 3
        }
      }
      cd$iter[i] <- iter
    }
    cd$iter[nrow(cd)] <- iter
    
    cdOut <- bind_rows(cdOut,cd)
    
    elapsed <- Sys.time() - startTime
    print(c(iter,elapsed))
  }  
  
  return(cdOut)
}
