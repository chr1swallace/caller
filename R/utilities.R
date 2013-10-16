beta.sd <- function(a,b) {                                        #  cat(a,b,"\n")
  sqrt(a*b/((a+b)^2 * (a+b+1)))
}


beta.mn <- function(a,b) {
  a/(a+b)
}
##' Given mean, sd of beta variable, calculate alpha, beta parameters
##'
##' @title getab.beta
##' @param mu mean
##' @param sigma sd
##' @return a list of alpha, beta vectors
##' @author Chris Wallace
getab.beta <- function(mu,sigma) {
  v <- sigma^2
  if(any(is.na(sigma)))
    v[is.na(sigma)] <- 0.01
  wh.low <- which(v==0 & mu<1e-16)
  if(length(wh.low)) {
    mu.in <- mu
    v.in <- v    
    mu <- mu[-wh.low]
    v <- v[-wh.low]
  }
  alpha=mu*(mu*(1-mu)/v - 1)
  beta=(1-mu)*(mu*(1-mu)/v - 1)

  ## deal with any negs
  wh.neg <- which(alpha<0 & beta<0)
  while(length(wh.neg)) {
    tmp <- getab.beta(mu[wh.neg],sqrt(v)[wh.neg]/10)
    alpha[wh.neg] <- tmp$alpha
    beta[wh.neg] <- tmp$beta
    wh.neg <- which(alpha<0 & beta<0)
  }
  
  ## add back in any low
  if(length(wh.low)) {
    a.tmp <- b.tmp <- numeric(length(mu.in))
    a.tmp[-wh.low] <- alpha
    b.tmp[-wh.low] <- beta
    a.tmp[wh.low] <- 1000*mu.in[wh.low]
    b.tmp[wh.low] <- 1000
    alpha <- a.tmp
    beta <- b.tmp
  }
  
  ##     cat("mean.obs=",mu,"mean.ab=",alpha/(beta+alpha),"\n")
  ##     cat("var.obs=",sigma^2,"var.ab=",alpha*beta/((alpha+beta)^2 * (alpha+beta+1)),"\n")
  list(alpha=alpha,beta=beta)
}
