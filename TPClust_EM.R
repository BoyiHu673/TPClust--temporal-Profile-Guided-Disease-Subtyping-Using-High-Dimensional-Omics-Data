for (k in 1:ncls){
  tmpw <- w[,k]
  vw <- rep(tmpw, DD)
  aa <- rep(alpha[,k], DD)
  aavar <- rep(alpha_var[,k], DD)
  rlambda <- tmp_rlam
  source('H2.R')
}
slambda <- tmp_slam
source('H1_sparse.R')
tmpu <- c()
for (k in 1:ncls){
  tmpu <- cbind(tmpu,uhat)
}
detach(package:CVXR,unload = TRUE)
library(CVXR)

old_that <- that
old_rhat <- rhat
old_ghat <- ghat
old_shat <- shat
ghat <- cbind(0,ghat)
tol <- 1e-4
est_diff <- 1
emcount <- 0
while (est_diff > tol & emcount<=60){
  w <- matrix(0,nrow=n,ncol=ncls)
  alpha <- matrix(0,nrow=n,ncol=ncls)
  alpha_var <- matrix(0,nrow=n,ncol=ncls)
  log_pi <- cbind(1,g)%*%ghat
  log_pi[log_pi>=700] <- 700
  pi <- exp(log_pi)/apply(exp(log_pi),1,sum)
  for (i in 1:n){
    tmp_log <- c()
    for (k in 1:ncls){
      mu <- 0
      bt <- getbasismatrix(t[i,!is.na(y[i,])&!is.na(t[i,])], basisobj)
      for (j in 1:J){
        mu <- mu + mtx[i,c(1:D)[!is.na(y[i,])&!is.na(t[i,])]+(j-1)*D]*bt%*%rhat[k,1:nb+(j-1)*nb]
      }
      for (l in 1:L){
        mu <- mu + datx[i,l]*bt%*%that[k,1:nb+(l-1)*nb]
      }
      sigma <- shat[k,2]^2*matrix(1,nrow=DD[i],ncol=DD[i]) + diag(shat[k,1]^2,DD[i])
      tmp_log <- c(tmp_log, log(dmvnorm(y[i,!is.na(y[i,])&!is.na(t[i,])],mean = mu, sigma = sigma)))
    }
    if (min(tmp_log)>0&sum(is.finite(tmp_log))==ncls) tmp_log <- tmp_log - min(tmp_log)
    for (k in 1:ncls){
      w[i,k] <- exp(tmp_log[k])*pi[i,k]
    }
  }
  w <- w/apply(w,1,sum)
  w[w<=1e-10] <- 1e-10
  for (i in 1:n){
    for (k in 1:ncls){
      mu <- 0
      alpha[i,k] <- 0
      alpha_var[i,k] <- 0
      bt <- getbasismatrix(t[i,!is.na(y[i,])&!is.na(t[i,])], basisobj)
      for (j in 1:J){
        mu <- mu + mtx[i,c(1:D)[!is.na(y[i,])&!is.na(t[i,])]+(j-1)*D]*bt%*%rhat[k,1:nb+(j-1)*nb]
      }
      for (l in 1:L){
        mu <- mu+ datx[i,l]*bt%*%that[k,1:nb+(l-1)*nb]
      }
      sigma <- shat[k,2]^2*matrix(1,nrow=DD[i],ncol=DD[i]) + diag(shat[k,1]^2,DD[i])
      alpha[i,k] <- w[i,k]*(shat[k,2]^2*rep(1,DD[i])%*%solve(sigma)%*%as.matrix(y[i,!is.na(y[i,])&!is.na(t[i,])] - mu))
      alpha_var[i,k] <- w[i,k]*(shat[k,2]^2 - shat[k,2]^4*rep(1,DD[i])%*%solve(sigma)%*%as.matrix(rep(1,DD[i])))
    }
  }
 
  for (k in 1:ncls){
    tmpw <- w[,k]
    vw <- rep(tmpw, DD)
    aa <- rep(alpha[,k], DD)
    aavar <- rep(alpha_var[,k], DD)
    rlambda <- tmp_rlam
    source('H2.R')
  }
  source('H1_sparse.R')
  detach(package:CVXR,unload = TRUE)
  library(CVXR)
  
  tmp_diff <- c(mean((that-old_that)^2), mean((rhat-old_rhat)^2), mean((ghat-old_ghat)^2), mean((shat-old_shat)^2))
  
  est_diff <- max(tmp_diff)
  old_that <- that
  old_rhat <- rhat
  old_ghat <- ghat
  old_shat <- shat
 
  ghat <- cbind(0,ghat)
  emcount <- emcount + 1
  if (emcount%%10==0){
    print(c('EM iterations:',emcount))
    print(c('Stopping criterion:', round(est_diff,6)))
    print(c('# of non-zero omics:', sum(abs(ghat)>=1e-3)))
  }
}
