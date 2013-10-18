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
##' @param max.lrr ignored, in future will allow maximization in theta only
##' @return a clusterfit object
##' @export
##' @author Chris Wallace
fit.em <- function(df.in,theta,lrr,tol=1e-2,verbose=TRUE,maxit=1e4, eps=1e-2,
                    use.deriv=TRUE, max.lrr=TRUE) {
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
  pars <- c(mu=df.in@R.mean[ df.in@R.opt ], sigma=df.in@R.sd[ df.in@R.opt ],
##             theta.mean=beta.mn(a=df.in@theta.a[ df.in@theta.opt ], b=df.in@theta.b[ df.in@theta.opt ]),
##             theta.sd=beta.sd(a=df.in@theta.a[ df.in@theta.opt ], b=df.in@theta.b[ df.in@theta.opt ]),
            a=df.in@theta.a[ df.in@theta.opt ], b=df.in@theta.b[ df.in@theta.opt ])
  
  ## define lots of functions within this function to use theta, lrr in this environment  
  
  parvec <- function(pars,index=TRUE) {
    a <- pars[grep("^a",names(pars))]
    b <- pars[grep("^b",names(pars))]
    mu <- pars[grep("^mu",names(pars))]
    sigma <- pars[grep("^sigma",names(pars))]
    names(a) <- sub(".*\\.","",names(a))
    names(b) <- sub(".*\\.","",names(b))
    names(mu) <- sub(".*\\.","",names(mu))
    names(sigma) <- sub(".*\\.","",names(sigma))
    if(index) {
      a <- a[as.character(dfsumm$theta.index)]
      b <- b[as.character(dfsumm$theta.index)]
      mu <- mu[as.character(dfsumm$R.index)]
      sigma <- sigma[as.character(dfsumm$R.index)]
      a[ dfsumm$copies==0 ] <- b[ dfsumm$copies==0 ] <- 1
      b[ dfsumm$copies==0 ] <- 1
      mu[ dfsumm$copies==0 ] <- dfsumm$R.mean[ dfsumm$copies==0 ]
      sigma[ dfsumm$copies==0 ] <- dfsumm$R.sd[ dfsumm$copies==0 ]
    }
    tm <- beta.mn(a,b)
    ts <- beta.sd(a,b)
    return(list(theta.mean=tm,theta.sd=ts,a=a,b=b,mu=mu,sigma=sigma))
  }

  pars.fail <- function(pars) {
  parv <- parvec(pars,index=FALSE)
  if(!identical(order(parv$theta.mean),seq_along(parv$theta.mean))) # preserve order of theta.mean
    return(TRUE)
  if("mu" %in% names(parv) && !identical(order(parv$mu),seq_along(parv$mu))) ## preserve order of LRR
    return(TRUE)
  if(any(parv$a < 0 | parv$b < 0)) # constraints imposed by beta distribution    
    return(TRUE)
  if(any(parv$theta.sd>0.1)) # keep theta clusters tight
    return(TRUE)
##   if(parv$theta.mean[2]<0.15 || parv$theta.mean[ length(parv$theta.mean)-1 ]>0.85) # middle clusters not too close to edge
##     return(TRUE)
  return(FALSE)
}
    
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
    ea <- eb <- matrix(0,length(theta),length(df.in@theta.a[ df.in@theta.opt ]),
                       dimnames=list(NULL,names(df.in@theta.a[ df.in@theta.opt ])))
    em <- es <- matrix(0,length(theta),length(df.in@R.mean[ df.in@R.opt ]),
                       dimnames=list(NULL,names(df.in@R.mean[ df.in@R.opt ])))
    for(i in 1:ngroup) {
      j <- as.character(dfsumm$theta.index[i])
      k <- as.character(dfsumm$R.index[i])
      tmp <- d.single(mu=parv$mu[i], sigma=parv$sigma[i], a=parv$a[i], b=parv$b[i], px[,i])
      if(dfsumm$theta.opt[i]) {
        ea[,j] <- ea[,j] + tmp$a
        eb[,j] <- ea[,j] + tmp$b
      }
      if(dfsumm$R.opt[i]) {
        em[,k] <- em[,k] + tmp$mu
        es[,k] <- es[,k] + tmp$sigma
      }
    }
    L <- lhood(pars, px, sumlog=FALSE)
    ret <- c(colSums(em/L),colSums(es/L),colSums(ea/L),colSums(eb/L))
    ret[ abs(ret)>1e100 ] <- sign(ret)[ abs(ret) > 1e100 ] * 1e100
    return(ret)
  }
   
  ## now start the work
  
  nit <- 0
  df <- 1
#  cols <- brewer.pal(ngroup, "Paired")
                                        #  hwe <- ngroup>1
  hwe <- FALSE
  value <- numeric(maxit)
  while(hwe | (df>tol & nit<maxit)) {
    
    nit <- nit+1    
    parv <- parvec(pars)
    
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
                      method="BFGS")
#      ,control=list(trace=verbose))
    } else {
      mx <- optim(par=pars,
                      fn=lhood,
                      px=px)
#      ,control=list(trace=verbose))
    }
    value[i] <- mx$value
    pars <- mx$par
    if(verbose) {
      message(nit,value[i],"\n")
      print(pars)
    }
  }
  subin <- function(old,opt,new) {
    ret <- old
    ret[ opt ] <- new
    return(ret)
  }
  df.out <- new("clusterdef",
                a1=df.in@a1,
                a2=df.in@a2,
                pi=colMeans(px),
                R.index=df.in@R.index,
                theta.index=df.in@theta.index,
                R.mean=subin(df.in@R.mean, df.in@R.opt, pars[ grep("^mu",names(pars)) ]),
                R.sd=subin(df.in@R.sd, df.in@R.opt, pars[ grep("^sigma",names(pars)) ]),
                theta.a=subin(df.in@theta.a, df.in@theta.opt, pars[ grep("^a",names(pars)) ]),
                theta.b=subin(df.in@theta.b, df.in@theta.opt, pars[ grep("^b",names(pars)) ]),
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

