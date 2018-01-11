# Burn injury data example
# simulations as in Huang (2016)
library(Matrix)

library(Rcpp)
sourceCpp("C:/Users/skm/Dropbox/mts/src/dvh.cpp")

burn_sim = function(N,show_plot=FALSE){
  b1 = c(6.698,0.0039)
  sig = 1.26*1.26
  gam = 5
  b2 = c(-4.0521,0.0527)
  ages = runif(N,20,80)
  X = cbind(1,ages)
  pis = 1/(1+exp(-X%*%b2))
  Y2 = rbinom(N,1,pis)
  mu1 = X%*%b1
  adj1 = gam*(Y2-pis)
  Y1 = rnorm(N,mu1+adj1,sqrt(sig))
  if(show_plot){
    yf = lm(Y1~adj1)
    R2 = sprintf("%.3f",cor(adj1,Y1)^2)
    plot(adj1,Y1,main=bquote(R^2==.(R2)))
    abline(a=yf$coef[1],b=yf$coef[2],col="steelblue4",lwd=2,lty=5) 
  }
  m1 = lm(Y1~ages)
  m2 = glm(Y2~ages,family=quasibinomial)
  m3 = glm(Y2~ages,family=binomial(link="logit"))
  Xmm = model.matrix(m1)
  s2 = sigma(m1)^2
  # multiplier bootstrap test
  bootTest = mbTest(cbind(X*as.numeric(Y1-X%*%b1),X*as.numeric(Y2-pis)),5000,c(-1,2))
  # GEE with working independence calculations
  muh1 = m1$fitted.values
  muh2 = m2$fitted.values
  res1 = Y1 - muh1
  res2 = Y2 - muh2
  W1 = rep(s2,N)
  W2 = summary(m2)$dispersion*muh2*(1-muh2)
  model_vcov = Matrix::bdiag(vcov(m1),vcov(m2))
  DW1 = sweep(Xmm,1,W1,`/`)
  DW2 = sweep(Xmm*(muh2*(1-muh2)),1,W2,`/`)
  DWS = cbind(sweep(DW1,1,res1,`*`),sweep(DW2,1,res2,`*`))
  meats = sapply(1:N,function(i){outer(DWS[i,],DWS[i,])},simplify = "array")
  meat = apply(meats,c(1,2),sum)
  # model-based sandwich estimator calculations
  naive_vcov = bdiag(vcov(m1),vcov(m3))
  sc1 = Xmm*(m1$residuals/s2)
  pih = 1/(1+exp(-Xmm%*%m3$coefficients))
  sc2 = Xmm*as.numeric(Y2-pih)
  SCZ = cbind(sc1,sc2)
  meats2 = sapply(1:N,function(i){outer(SCZ[i,],SCZ[i,])},simplify = "array")
  meat2 = apply(meats2,c(1,2),sum)
  list("model" = 
         list("beta" = c(m1$coefficients,m3$coefficients),
              "sandwich" = naive_vcov %*% meat2 %*% naive_vcov,
              "naive" = naive_vcov),
       "gee" =
         list("beta" = c(m1$coefficients,m2$coefficients),
              "sandwich" = model_vcov %*% meat %*% model_vcov,
              "naive" = model_vcov),
       "mbt.p" = bootTest$p.value)
}

# example:
(BS = burn_sim(200))

# the F test for a given contrast matrix M
# bcov is (typically) the sandwich estimator
f_test = function(beta,bcov,delta,M,N){
  r = nrow(M)
  stopifnot(length(beta)==nrow(bcov) && nrow(bcov)==ncol(bcov))
  stopifnot(length(delta) == r && ncol(M)==nrow(bcov))
  t1 = M%*%beta - delta
  t2 = M%*%bcov%*%t(M)
  Fstat = as.numeric((N/r)*t(t1)%*%solve(N*t2)%*%t1) # cast from Matrix S4 member
  list("F" = Fstat,"p.value" = 1-pf(Fstat,r,N-length(beta)))
}

dd1 = c(6.698,0.0039,-4.0521,0.0527)
M1 = diag(4) # for testing everything against true values ("type I error")
M2 = rbind(c(0,1,0,0),c(0,0,0,1)) # testing 'joint effect' of slopes
MC = 10000  # number of Monte Carlo simulations per sample size
NN = c(50,100,200,400,800) # sample sizes

for(N in NN){
  cat(paste0("\n\n\n\nN = ",N,"\n\n"))
  FITS1 = vector("list",MC)
  FITS2 = vector("list",MC)
  FITS3 = matrix(1,nrow=MC,ncol=2)
  # F-stat and p value for both sandwich and naive estimators
  Fstats1 = matrix(0,nrow=MC,ncol=4)
  Fstats2 = matrix(0,nrow=MC,ncol=4)
  i = 1
  while(i <= MC){
    tryCatch({
      fit = burn_sim(N,FALSE)
      fG = fit$gee; fM = fit$model
      FITS1[[i]] = fG; FITS2[[i]] = fM; FITS3[i,] = fit$mbt.p
      Fstats1[i,1:2] = unlist(f_test(fG$beta,fG$sandwich,dd1,M1,N))
      Fstats1[i,3:4] = unlist(f_test(fG$beta,fG$naive,dd1,M1,N))
      Fstats2[i,1:2] = unlist(f_test(fM$beta,fM$sandwich,dd1,M1,N))
      Fstats2[i,3:4] = unlist(f_test(fM$beta,fM$naive,dd1,M1,N))
      # cat(paste0(i,'\n'))
      NULL
      },
      error = function(e){
        return(NA)
    })
    i <<- i+1
    if(i%%500 == 0) cat(paste0(i,"..."))
  }
  
  cat('\n')
  print("GEE type I errors:")
  print(sapply(c(0.01,0.05,0.10),function(alpha){mean(Fstats1[,2] < alpha)}))
  print("Model type I errors:")
  print(sapply(c(0.01,0.05,0.10),function(alpha){mean(Fstats2[,2] < alpha)}))
  print("mbt type I errors:")
  print(sapply(c(0.01,0.05,0.10),function(alpha){colMeans(FITS3 < alpha)}))
  cor1 = sapply(FITS1,function(e){cov2cor(e$sandwich)[4,2]})
  cor2 = sapply(FITS2,function(e){cov2cor(e$sandwich)[4,2]})
  print("Correlations:")
  print(paste0(sprintf("%.3f",mean(cor1))," (",sprintf("%.3f",sd(cor1)),")"))
  print(paste0(sprintf("%.3f",mean(cor2))," (",sprintf("%.3f",sd(cor2)),")"))
}



