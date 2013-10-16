setMethod("summary","clusterdef",
          definition=function(object) {
            sd <- beta.sd(object@theta.a,object@theta.b)
            mn <- beta.mn(object@theta.a,object@theta.b)
            data.frame(a1=object@a1,
                             a2=object@a2,
                             copies=object@a1 + object@a2,
                            pi=object@pi,
                             R.index=object@R.index,
                             R.mean=object@R.mean[ object@R.index ],
                             R.sd=object@R.sd[ object@R.index ],
                       R.opt=object@R.opt[ object@R.index ],
                             theta.index=object@theta.index,
                             theta.mean=mn[ object@theta.index ],
                             theta.sd=sd[ object@theta.index ],
                             theta.a=object@theta.a[ object@theta.index ],
                             theta.b=object@theta.b[ object@theta.index ],
                       theta.opt=object@theta.opt[ object@theta.index ]
                            )
          })

setMethod("show","clusterdef",
          definition=function(object) {
            cat("clusterdef object defining",length(object@a1),"clusters with 2 x",length(object@R.mean),"R parameters and 2 x",length(object@theta.a),"theta parameters.\n")
            show(summary(object)[,c("a1","a2","copies","pi","R.mean","R.sd","theta.mean","theta.sd")])
          })
                             
                             
setMethod("show","clusterfit",
          definition=function(object) {
            show(object@clusters)
            cat("Fitted to",length(object@theta),"observations in",length(object@fit.info),"steps\n")
          })
                             
                             
          
              
    
              