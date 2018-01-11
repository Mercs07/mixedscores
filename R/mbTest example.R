# Example running the multiplier bootstrap test
# visit the README file for documentation

library(Rcpp)
sourceCpp("https://github.com/smarches/mixedscores/blob/master/src/dvh.cpp") # or your local path

z = cbind(rnorm(100,0.25,1),rnorm(100,0,2),rexp(100,1)-1)

res = mbTest(z,1000,c(-1,2,1))

res$p.value
res$obs.norm
head(res$norms)
matplot(res$norms,type="p")
cor(res$norms)