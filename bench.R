#source("newallelic.R")
library("newallelic")

systematicTestNewAllelic <- function(n1,n2,incr) {
  print( c("n1",n1,"n2",n2,"incr",incr)  );
  table <- matrix(c(0,0,0,0,0,0),nrow=2,ncol=3, byrow=TRUE)
  pvalues <- c()
  fisher22 <- c()
  fisher23 <- c()
  table22 <- matrix(0,nrow=2,ncol=2, byrow=TRUE)
  for(a in seq(0,n1,by=incr)) {
    print(c("a=",a))
    table[1,1] <- a
    for(b in seq(0,n1-a,by=incr)) {
      table[1,2] <- b
      table22[1,1] <- a + a + b
      table[1,3] <- n1-a-b
      table22[1,2] <- table[1,3]*2 + b
      for(d in seq(0,n2,by=incr)) {
        table[2,1] <- d
       
        for(e in seq(0,n2-d,by=incr)) {
          table[2,2] <- e
          table[2,3] <- n2-d-e
          table22[2,1] <- d + d + e
          table22[2,2] <- table[2,3]*2 + e
#          fisher<-Fisher22(table22,hybrid = FALSE,conf.int = FALSE)$p.value
#          f <- Fisher23(table,hybrid = FALSE,conf.int = FALSE)$p.value
#          f <- fisher.test(table,hybrid = FALSE,conf.int = FALSE)$p.value
#          fisher22[length(fisher22)+1] <- fisher
#          fisher23[length(fisher23)+1] <- f
          p <- newallelic_nocheck(table)
   #       print( p )
          pvalues[length(pvalues)+1] <- p
        }
      }
    }
  }
  
  print( c("Got ",length(pvalues)," pvalues"))
  print( c("Sum=",sum(pvalues)))

  return( list(fisher22=fisher22,newallelic=pvalues,fisher23=fisher23))
}

#Rprof("newallelic.out")
l <- systematicTestNewAllelic(300,350,80)
#Rprof(NULL)
print(l$pvalues)
