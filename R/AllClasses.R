#' @exportClass clusterdef
setClass("clusterdef",
         slots=c(a1="numeric",
           a2="numeric",
           pi="numeric",
           R.index="numeric",
           theta.index="numeric",
           R.mean="numeric",
           R.sd="numeric",
           theta.a="numeric",
           theta.b="numeric",
           R.opt="logical",
           theta.opt="logical"),
         validity=function(object) {
           n <- length(object@a1)
           if(length(object@a2)!=n ||
              length(object@R.index)!=n ||
              length(object@theta.index)!=n ||
              length(object@pi)!=n)
             stop("a1, a2, R.index, theta.index, pi need to be of equal length")
           m <- max(object@R.index, na.rm=TRUE)
           if(length(object@R.mean)!=m ||
              length(object@R.sd)!=m)
             stop("not all indices in R.index given in R.mean, R.sd")
           m <- max(object@theta.index, na.rm=TRUE)
          if(length(object@theta.a)!=m ||
              length(object@theta.b)!=m)
             stop("not all indices in theta.index given in theta.a, theta.b")
         })


#' @exportClass clusterfit
setClass("clusterfit",
         slots=c(theta="numeric",
             lrr="numeric",
             pp="matrix",
             fit.info="numeric",
           clusters="clusterdef"),
         validity=function(object) {
           n <- length(object@theta)
           if(length(object@lrr)!=n ||
              nrow(object@pp) != n)
             stop("theta and lrr must have length == nrow(pp)")
         })


