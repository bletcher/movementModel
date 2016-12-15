
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