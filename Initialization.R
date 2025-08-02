# initial weight matrix
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