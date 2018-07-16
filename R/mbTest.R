# R only version of mbTest

lp_norm = function(x,i){
  if(i==-1) max(abs(x)) else (sum(abs(x)^i))^(1/i)
}

mb_test = function(x,B,norms){
  norms = sort(unique(as.integer(norms)))
  norms = norms[norms == -1 || norms > 0]
  V = sqrt(colSums(x*x))
  Bnorms = sapply(1:B,function(i){
    Utilde = (t(x)%*%(2*rbinom(NROW(x),1,0.5)-1))/V
    sapply(norms,lp_norm,"x"=Utilde)
  })
  Onorm = sapply(norms,lp_norm,"x"= colSums(x)/V,simplify="array")
  pv = if(length(norms)==1) mean(Onorm <= Bnorms) else apply(Onorm <= Bnorms,1,mean)
  list("obs.norm" = Onorm,
       "p.value" = pv,
       "norms" = t(Bnorms))
}

# library(MASS)
# Y = mvrnorm(100,c(1,-1),matrix(c(1,-0.7,-0.7,1),ncol=2))
# plot(mbTest(Y,2000,c(-1,2))$norms)
# plot(mbTest(cbind(runif(100),runif(100)),2000,c(-1,2))$norms)
# hist(mbTest(Y,2000,c(-1,2))$norms[,1])
# hist(mbTest(cbind(runif(100),runif(100)),2000,c(-1,2))$norms[,1],col="salmon",add=TRUE)
