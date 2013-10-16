if(FALSE) {
  test.data <- list(Z=c(0.1,0.1,0.1,0.4,0.4,0.6, 0.8,0.7,0.6),s=c(10,3,20,4,50,2,50,40,10))
  source("R/mix-betabinomial.R")
  bb.result <- mixbb.fit(test.data$Z,test.data$s,ngroup=3)
  beta.result <- mixbb.fit(test.data$Z,test.data$s,ngroup=3,bb=FALSE)
  bb.result$par
  beta.result$par
}


mixbb.fit <- function(Z,s,ngroup=3,bb=TRUE,...) {
  ## overall fitting function
  ## fit a beta/betabinomial by EM
  ## Z is the proportion, s the number of observations per sample, ngroup the number of mixtures
  ## bb=TRUE fits a betabinomial, which makes use of s
  ## bb=FALSE fits a mix of betas to Z only
  
  pars <- mix.getpars(Z,ngroup=ngroup)

  ## get alpha,beta for each group - limit as n-> infty
  alphbet <- beta.ab(pars$mu.est,pars$sigma.est)
  alph <- alphbet$alpha
  bet <- alphbet$beta
  p <- c(pars$pi,1-sum(pars$pi))

  if(ngroup==1) {
    theta <- c(a1=alph,b1=bet)
  } else {
    theta <- unlist(list(p=p[1:(ngroup-1)],a=alph,b=bet))
  }

  mix <- mixbb.em(theta,Z,s,bb=bb,...)

  mnsd <- beta.mnsd(mix$par[grep("a",names(mix$par))],
                    mix$par[grep("b",names(mix$par))])
  mix$group.mean <- mnsd$mn
  mix$group.sd <- mnsd$sd
  mix$posterior.prob <- mixbb.call(Z,s,mix$par,bb=bb)
  mix$posterior.best <- apply(mix$posterior.prob,1,max)
  mix$group.best <- apply(mix$posterior.prob,1,which.max)
  names(mix$group.best) <- names(mix$prob.best) <- rownames(mix$posterior.prob) <- names(s)
  
  return(mix)
}

################################################################################

mixbb.call <- function(Z,s,theta,bb=TRUE) {
  p <- theta[grep("p",names(theta))]
  if(!length(p)) {
    p <- 1
  }
   ## alpha, beta from the underlying beta distribution
  alph <- theta[grep("a",names(theta))]
  bet <- theta[grep("b",names(theta))]
  ngroup <- length(alph)
  
  ## generate probs to make calls
  px <- matrix(1,length(Z),ngroup)
  
  if(bb) {  ## alpha, beta from beta approx of beta binomial
    S <- ((alph+bet+1) %o% s)/(outer(alph+bet,s,"+")) - 1
    a <- S*(alph/(alph+bet))
    b <- S*(bet/(alph+bet))
  }
  
  for(i in 1:ngroup) {
    if(bb)
      px[,i] <- beta.lhood(Z,a[i,],b[i,],p[i])
    else
      px[,i] <- beta.lhood(Z,alph[i],bet[i],p[i])
  }
  px <- px/rowSums(px)
  
  return(px)
}

################################################################################

mixbb.em <- function(theta.in,Z,s,## limits=NA,
                    tol=1e-4,verbose=FALSE,maxit=1e4,bb=TRUE) {
  ## fit a beta/betabinomial by EM
  ## Z is the proportion, s the number of observations per sample
  ## theta is a list of alpha, beta and p for the currently fitted mixture
  ## finds the MLE of alpha, beta given fixed p
  mx <- NULL
  p <- theta.in[grep("p",names(theta.in))]
  if(length(p)) {
    p <- c(p,1-sum(p))
  } else {
    p <- 1
  }
  ngroup <- length(p)
  
  ## alpha, beta from the underlying beta distribution
  alph <- theta.in[grep("a",names(theta.in))]
  bet <- theta.in[grep("b",names(theta.in))]
    
  ## P(group[i]=j | Z[i],n[i]) = P(Z[i]|group[i]==j,n[i]) * P(group==j) / P(Z[i]|n[i])
  px <- matrix(1,length(Z),ngroup)
  
  nit <- 0
  df <- 1

  while(df>tol & nit<maxit) {
    
    nit <- nit+1
    
    if(bb) { ## alpha, beta from beta approx of beta binomial
      S <- ((alph+bet+1) %o% s)/(outer(alph+bet,s,"+")) - 1
      a <- S*(alph/(alph+bet))
      b <- S*(bet/(alph+bet))
    }  # else just do standard beta

    ## E step
    if(verbose) { cat("E step\n") }
    if(ngroup>1) { # px=1 for ngroup==1
      px.old <- px
      for(i in 1:ngroup) {
        if(bb)
          px[,i] <- beta.lhood(Z,a[i,],b[i,],p[i])
        else
          px[,i] <- beta.lhood(Z,alph[i],bet[i],p[i])
      }
      px <- px/rowSums(px)
      if(any(is.nan(px)))
        px[is.nan(px)] <- 0
      df <- max(abs(px-px.old))
    } else {
      df <- 0 # no point iterating over an E step that doesn't change
    }
    
    ## M step
    if(verbose) { cat("M step\n") }
    p <- colMeans(px)
    mx <- try(optim(par=c(alph,bet),
                    fn=mixbb.lhood,
                    p=px,Z=Z,s=s,## limits=limits,
                    bb=bb))
    
    if(inherits(mx,"try-error")) {
      mx <- NULL
      next
    }
    
    ## alpha, beta from the underlying beta distribution
    alph <-  mx$par[grep("a",names(mx$par))]
    bet <- mx$par[grep("b",names(mx$par))]
  }

  if(is.null(mx))
    return(mx)
  
  if(ngroup==1) {
    mx$par <- c(mx$par,p1=1)
  } else {
    mx$par <- c(mx$par,p=unlist(p))
  }
  return(mx)
}

################################################################################


beta.lhood <- function(Z,a,b,p) {
  tmp <- try(dbeta(Z,shape1=a,shape2=b) * p)
  if(inherits(tmp,"error")) {
    tmp <- rep(1e-64,length(Z)) # not 0, as we may want to log this later
  }
  if(length(wh <- which(tmp==0))) { # avoid numerical problems later
    tmp[wh] <- 1e-64
  }
  if(length(wh <- which(is.infinite(tmp)))) { # avoid numerical problems later
    tmp[wh] <- 1e64
  }
  return(tmp)
}

mixbb.lhood <- function(theta,p,Z,s,## limits=NA,
                          bb=TRUE) {
  ## calulates lhood for either a beta or betabinomial (if bb is TRUE) mixture distribution
  ## Z is the proportion, s the number of observations per sample
  ## theta is a list of alpha, beta and p for the currently fitted mixture
 
  if(any(theta<0))
    return(NA)

  ## alpha, beta from the underlying beta distribution
  alph <- theta[grep("a",names(theta))]
  bet <- theta[grep("b",names(theta))]

  if(bb) { ## alpha, beta from beta approx of beta binomial
    S <- ((alph+bet+1) %o% s)/(outer(alph+bet,s,"+")) - 1
    a <- S*(alph/(alph+bet))
    b <- S*(bet/(alph+bet))
  }  # else just do standard beta

  ngroup <- length(alph)
  e <- rep(0,length(Z))
  for(i in 1:ngroup) {
    if(bb)
      e <- e+beta.lhood(Z,a[i,],b[i,],p[,i])
    else
      e <- e+beta.lhood(Z,alph[i],bet[i],p[,i])
  }
  ## means should be ordered, and central group should be in the middle
  mn <- beta.mn(alph,bet)
    if(ngroup>=2 && !identical(mn,sort(mn)))
    return(NA)
  sd <- beta.sd(alph,bet)
  -2 * sum(log(e))
}

################################################################################

## helper functions to convert alpha, beta of beta distributions to mean, sd and vice versa

beta.sd <- function(a,b) {  sqrt(a*b/((a+b)^2 * (a+b+1))) }
beta.mn <- function(a,b) {  a/(a+b) }
beta.mnsd <- function(a,b,par=NULL) {
  if(!is.null(par)) {
    a <- par[grep("a",names(par))]
    b <- par[grep("b",names(par))]
  }
  list(mn=beta.mn(a,b),sd=beta.sd(a,b))
}
beta.ab <- function(mu,sigma) {
  ## returns alpha, beta parameters for a beta distribution with mean mu, sd sigma
  v <- sigma^2
  ## can reach some strange values of sigma.  drag these back to within normal limits.
  if(length(wh <- which(is.na(sigma) | sigma==0)))
    v[wh] <- 0.01
  alpha=mu*(mu*(1-mu)/v - 1)
  beta=(1-mu)*(mu*(1-mu)/v - 1)

  ## deal with any negs by decreasing sigma and trying again
  wh.neg <- which(alpha<0 & beta<0)
  while(length(wh.neg)) {
    tmp <- beta.ab(mu[wh.neg],sqrt(v)[wh.neg]/10)
    alpha[wh.neg] <- tmp$alpha
    beta[wh.neg] <- tmp$beta
    wh.neg <- which(alpha<0 & beta<0)
  }

  list(alpha=alpha,beta=beta)
}

