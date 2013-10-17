##' Fit a beta mixture distribution to theta, assuming LRR follows an established Gaussian mixture distribution
##'
##' @title fit.em
##' @param df.in clusterdef object
##' @param theta theta
##' @param lrr LRR
##' @param tol how small a change in group membership probabilities prompts continued optimization
##' @param verbose print messages if TRUE
##' @param maxit maximum number of iterations
##' @param eps theta cannot be exactly 0 or 1, and values less than eps from 0 or 1 are set to eos and 1-eps respectively
##' @param use.deriv if TRUE, use gradients to optimize quicker.  Currently off by default because it DOESN'T WORK!
##' @return a clusterfit object
##' @export
##' @author Chris Wallace
both.em <- function(df.in,theta,lrr,tol=1e-2,verbose=TRUE,maxit=1e4, eps=1e-2, use.deriv=TRUE, max.lrr=FALSE) {
  mx <- NULL
  p <- df.in@pi
  ngroup <- length(p)
  
  ## theta \in (0,1)
  theta[ theta<=eps ] <- eps
  theta[ theta>=1-eps ] <- 1-eps
  
  dfsumm <- summary(df.in)
  
  ## probabilities of group membership
  px <- matrix(dfsumm$pi,length(theta),ngroup,byrow=TRUE)
  ## parameter vector
  pars <- c(mu=df.in@R.mean, sigma=df.in@R.sd,
            a=df.in@theta.a[ df.in@theta.opt ], b=df.in@theta.b[ df.in@theta.opt ])
  
  ## define lots of functions within this function to use theta, lrr in this environment  
  
  ## likelihood for a single group
  lhood.single <- function(mu, sigma, a, b, pi) {
    exp(dbeta(theta,shape1=a,shape2=b,log=TRUE) + dnorm(lrr,mean=mu,sd=sigma,log=TRUE) + log(pi))
  }
  
  
  ## likelihood function to be maximized
  lhood <- function(pars, px, sumlog=TRUE) {
    if(pars.fail(pars))
      return(NA)
    parv <- parvec(pars)
    ngroup <- length(parv$a)
    e <- numeric(length(theta))
    for(i in 1:ngroup) {
      ##    cat(i, (dbeta(theta,a[i],b[i]) * dnorm(lrr,mu[i],sigma[i]) * p[,i])[wh], "\n")
      e <- e + lhood.single(mu=parv$mu[i], sigma=parv$sigma[i], a=parv$a[i], b=parv$b[i], pi=px[,i])
    }
    if(!any(is.na(e)) & any(e==0)) {
      wh <- which(e==0)
      e[wh] <- 1e-64
    }
    if(!any(is.na(e)) & any(is.infinite(e))) {
      wh <- which(is.infinite(e))
      e[wh] <- 1e64
    }
    if(sumlog) {
      return(-sum(log(e)))
    } else {
      return(-e)
    }
  }
  
  d.single <- function(mu, sigma, a, b, pi, inf.value=1e400) {
#    B <- lbeta(a,b)
    C <- dnorm(lrr,mean=mu,sd=sigma,log=TRUE) + dbeta(theta, a, b, log=TRUE) + log(pi)
#      log(pi) + (b - 1)*log(1-theta) + (a-1)*log(theta) - B
    tmp <- list(a= exp(C) * (log(theta) - digamma(a) + digamma(a+b)),
                b= exp(C) * (log(1-theta) - digamma(b) + digamma(a+b)),
                mu=exp(C) * (lrr-mu)/sigma^2,
                sigma=exp(C) * ((lrr-mu)^2/sigma^2 - 1))
    return(tmp)
  }
  
  deriv <- function(pars, px) {
    parv <- parvec(pars)    
    ngroup <- length(parv$a)
    deriv.a <- deriv.b <- numeric(length(pars)/2)
    ea <- eb <- matrix(0,length(theta),length(grep("^a",names(pars))))
    em <- es <- matrix(0,length(theta),length(grep("mu",names(pars))))
    for(i in 1:ngroup) {
      j <- dfsumm$theta.index[i]
      k <- dfsumm$R.index[i]
      tmp <- d.single(mu=parv$mu[i], sigma=parv$sigma[i], a=parv$a[i], b=parv$b[i], px[,i])
      if(dfsumm$theta.opt[i]) {
        ea[,j] <- ea[,j] + tmp$a
        eb[,j] <- ea[,j] + tmp$b
      }
      em[,k] <- em[,k] + tmp$mu
      es[,k] <- es[,k] + tmp$sigma
    }
    L <- lhood(pars, px, sumlog=FALSE)
    ret <- c(colSums(em/L),colSums(es/L),colSums(ea/L),colSums(eb/L))
    ret[ abs(ret)>1e100 ] <- sign(ret)[ abs(ret) > 1e100 ] * 1e100
    return(ret)
  }
  
  ## now start the work
  
  nit <- 0
  df <- 1
  cols <- brewer.pal(ngroup, "Paired")
                                        #  hwe <- ngroup>1
  hwe <- FALSE
  value <- numeric(maxit)
  while(hwe | (df>tol & nit<maxit)) {
    
    nit <- nit+1    
    ab <- abvec(pars)    
    ms <- msvec(pars)    
    
    ## E step
    if(ngroup>1) { # px=1 for ngroup==1
      px.old <- px
      p <- colMeans(px)
      for(i in 1:ngroup) {
        px[,i] <- lhood.single(mu=parv$mu[i], sigma=parv$sigma[i], a=parv$a[i], b=parv$b[i], pi=p[i])
      }
      px <- px/rowSums(px) ## normalise
      if(any(is.nan(px)))
        px[is.nan(px)] <- 0
      df <- max(abs(px-px.old))
    } else {
      df <- 0 # no point iterating over an E step that doesn't change
    }
    
    ## M step
    p <- colMeans(px)
   if(use.deriv) {
      mx <- optim(par=pars,
                      fn=lhood,
                      gr=deriv,
                      px=px,
                      method="BFGS",
                  control=list(trace=verbose))
    } else {
      mx <- optim(par=pars,
                      fn=lhood,
                      px=px,
                      control=list(trace=verbose))
    }
    value[i] <- mx$value
    pars <- mx$par
    if(verbose) {
      cat(nit,value[i],"\n")
      print(pars)
    }
  }
  newa <- df.in@theta.a
  newa[ df.in@theta.opt ] <- pars[ grep("a",names(pars)) ]
  newb <- df.in@theta.b
  newb[ df.in@theta.opt ] <- pars[ grep("b",names(pars)) ]  
  df.out <- new("clusterdef",
                a1=df.in@a1,
                a2=df.in@a2,
                pi=colMeans(px),
                R.index=df.in@R.index,
                theta.index=df.in@theta.index,
                R.mean=df.in@R.mean,
                R.sd=df.in@R.sd,
                theta.a=newa,
                theta.b=newb,
                R.opt=df.in@R.opt,
                theta.opt=df.in@theta.opt)
  ##   df.out <- df.in
  ##   df.out@theta.a[ df.out@theta.opt ] <- mx$par[ grep("a",names(mx$par)) ]
  ##   df.out@theta.b[ df.out@theta.opt ] <- mx$par[ grep("b",names(mx$par)) ]
  ##   df.out@pi <- colMeans(px)
  return(new("clusterfit",
             theta=theta,
             lrr=lrr,
             pp=px,
             fit.info=value[1:nit],
             clusters=df.out))
}

