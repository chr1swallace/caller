##' Create a clusterdef object, specifying only the maximum number of total copies observed
##'
##' Default values for LRR means and sd are taken from ???.
##' The LRR_SD is quite dependent on data quality and normalisation but the
##' LRR-Mean should be fairly good.
##' 
##' LRR_mean:
##' Copy 0        1        2        3       4+
##' -3.527211 -0.664184 0.000000 0.395621 0.678345
##' LRR_sd:
##' Copy 0        1        2        3       4+
##' 1.329152 0.284338 0.159645 0.209089 0.191579                                 
##' @title clusterdef.ncopies
##' @param n maximum number of copies
##' @param R.mean default values of LRR means for 0..4 copies
##' @param R.sd default values of LRR sds for 0..4 copies
##' @param eps distance from 0 or 1 to assign theta values with theoretical values of 0 or 1
##' @return a clusterdef object
##' @export
##' @author Chris Wallace
clusterdef.ncopies <- function(n,
                               R.mean=c("0"=-3.527211,"1"=-0.664184,"2"=0.000000,"3"=0.395621,"4"=0.678345),
                               R.sd=c("0"=0.5 ,"1"=0.284338,"2"= 0.159645,"3"= 0.209089,"4"= 0.191579),
                               # NB reduced first R.sd from 1.329152 -> 0.5
                               eps=1e-2) {
## find all combinations of a1, a2 st 0 <= ai <= n and a1+a2 <= n
  a1 <- a2 <- 0:n
  df <- expand.grid(a1=a1,a2=a2)
  df <- df[ rowSums(df)<=n, ]
  ncopies <- df$a1+df$a2
  tmn <- atan2(df$a1,df$a2)/(pi/2)
  tmn[ df$a1==0 & df$a2==0 ] <- 0.5
  tmn[ tmn==0 ] <- eps
  tmn[ tmn==1 ] <- 1-eps
  tsd <- rep(0.01,nrow(df))
  tsd[ df$a1==0 & df$a2==0 ] <- sqrt(1/12)
  ab <- getab.beta(tmn,tsd)

  ## index theta parameters, setting NA for 0 copies
  ti <- df$a1/df$a2
  ## sort ab so that theta.index corresponds to correct value
  theta.index <- as.numeric(as.factor(ti))
  wh <- which(!duplicated(theta.index))
  a <- structure(ab$alpha[wh][ order(theta.index[wh]) ],names=sort(theta.index[wh]))
  b <- structure(ab$beta[wh][ order(theta.index[wh]) ],names=sort(theta.index[wh]))
  topt <- ifelse(df$a1==0 & df$a2==0, FALSE, TRUE)[wh][ order(theta.index[wh]) ]

  ## index R parameters, setting NA for 0 copies
  R.index <- pmin(ncopies,4) + 1
#  R.index[ ncopies==0 ] <- max(R.index) + 1
  wh <- which(!duplicated(R.index))
  m <- structure(R.mean[R.index[wh]][ order(R.index[wh]) ],names=sort(R.index[wh]))
  s <- structure(R.sd[R.index[wh]][ order(R.index[wh]) ],names=sort(R.index[wh]))
  ropt <- ifelse(df$a1==0 & df$a2==0, FALSE, TRUE)[wh][ order(R.index[wh]) ]
  
  new("clusterdef",
      a1=df$a1,
      a2=df$a2,
      pi=rep(1/nrow(df), nrow(df)),
      R.index=R.index,
      theta.index=theta.index,
      R.mean=m,
      R.sd=s,
      theta.a=a,
      theta.b=b,
      R.opt=ropt,
      theta.opt=topt)
}


initialize <- function(clusters,theta,lrr) {

  tmn <- beta.mn(clusters@theta.a[ clusters@theta.opt ], clusters@theta.b[ clusters@theta.opt ])
  tsd <- beta.sd(clusters@theta.a[ clusters@theta.opt ], clusters@theta.b[ clusters@theta.opt ])

  ## recenter middle clusters
  nmid <- length(clusters@theta.a[ clusters@theta.opt ]) - 2
  km <- kmeans(theta[ theta>0.15 & theta<0.85 ], nmid)
  wh <- 1 + (1:nmid)
  tmn[ wh ] <- sort(km$centers)
  tsd[ wh ] <- tapply(theta[ theta>0.15 & theta<0.85 ],km$cluster,sd)[order(km$centers)]

  ## add outside clusters
  tmn[ 1 ] <- mean(theta [ theta < 0.1 ])
  tsd[ 1 ] <- sd(theta [ theta < 0.1 ])
  n <- length(tmn)
  tmn[ n ] <- mean(theta [ theta > 0.9 ])
  tsd[ n ] <- sd(theta [ theta > 0.9 ])
  
  ## put it back in clusters
  ab <- getab.beta(tmn,tsd)
  clusters@theta.a[ clusters@theta.opt ] <- ab$alpha
  clusters@theta.b[ clusters@theta.opt ] <- ab$beta

  clusters
  
}
