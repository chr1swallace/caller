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

setMethod("summary","clusterfit",
          definition=function(object) {
            best.group <- apply(object@pp,1,which.max)
            best.prob <- apply(object@pp,1,max)
            data.frame(theta=object@theta,LRR=object@lrr,best.group=best.group,prob=best.prob)
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
                             
                             
setMethod("plot","clusterdef",
          definition=function(x,y) {
            dfsumm <- summary(x)
            dfsumm <- within(dfsumm, {
              theta.low <- qbeta(0.025, theta.a, theta.b)
              theta.high <- qbeta(0.975, theta.a, theta.b)
              R.low <- qnorm(0.025,R.mean,R.sd)
              R.high <- qnorm(0.975,R.mean,R.sd)              
            })
            ggplot(dfsumm,aes(col=as.factor(1:nrow(dfsumm)))) +
              geom_point(aes(x=theta.mean,y=R.mean)) +
                geom_segment(aes(x=theta.low,xend=theta.high,y=R.mean,yend=R.mean),alpha=0.5,
                             arrow=arrow(length = unit(0.1,"cm"), angle=90, ends="both")) +
                geom_segment(aes(x=theta.mean,xend=theta.mean,y=R.low,yend=R.high),alpha=0.5,lineend="round",
                             arrow=arrow(length = unit(0.1,"cm"), angle=90, ends="both"))
          })
              
setMethod("plot","clusterfit",
          definition=function(x,y) {
            groups <- summary(x)
            plot(x@clusters) + geom_point(data=groups, aes(x=theta,y=LRR,col=as.factor(best.group)),size=0.1)
          })
              
  
              
