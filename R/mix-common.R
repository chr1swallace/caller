mix.getpars <- function(Z,ngroup=3,inits=NULL,nstart=25) {
  ## get sensible starting values for a beta mixture of Z
  ## they will be guessed from nstart initiations of kmeans
  ## optionally with starting values for the means if supplied in
  ## inits$means which should be a numerical vector of length ngroup

  Z <- Z[!is.nan(Z) & !is.na(Z)]
  eps <- 1e-4
  
  if(ngroup>1) {
    if(!is.null(inits)) {
      km <- try(kmeans(Z,centers=inits$mean[1:ngroup],nstart=nstart))
    }
    if(is.null(inits) || inherits(km,"try-error")) {
      km <- try(kmeans(Z,centers=ngroup,nstart=nstart))
    }
    means.est <- km$centers[,1]
    sigma.est <- sapply(1:ngroup,function(i) sd(Z[km$cluster==i],na.rm=TRUE))
    sigma.est[is.na(sigma.est)] <- 0.01 # single obs in a group
    sigma.est[sigma.est==0] <- 0.01 # all obs equal
    pi <- km$size/length(Z)

    ord <- order(means.est)
    means.est <- means.est[ord]
    sigma.est <- sigma.est[ord]
    pi <- pi[ord]
  } else {
    means.est <- mean(Z,na.rm=TRUE)
    sigma.est <- sd(Z,na.rm=TRUE)
    pi <- 1
  }
  names(means.est) <- NULL
  
  return(pars=list(mu=if(is.null(inits)) { NA } else { inits$means },
           sigma=if(is.null(inits)) { NA } else { inits$sigma },
           mu.est=means.est,
           sigma.est=sigma.est,
           pi=pi))
}

################################################################################
