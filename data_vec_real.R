vy <- c(t(y))
ty <- t+y
ty <- c(t(ty))
vg <- matrix(0,nrow=length(vy),ncol=ncol(g))
for (j in 1:ncol(g)){
  vg[,j] <- rep(g[,j],each = D)
}
vg <- vg[is.na(vy)==FALSE,]
tmpxx <- c()
for (l in 1:L){
  tmpx <- c()
  for (i in 1:n){
    bt <- getbasismatrix(t[i,!is.na(y[i,])&!is.na(t[i,])], basisobj)
    tmpx <-rbind(tmpx, datx[i,l]*bt)
  }
  tmpxx <- cbind(tmpxx,tmpx)
}
xx <- tmpxx

if (!is.na(J)){
  mtxx <- c()
  for (j in 1:J){
    tmp_mtxx <- c()
    for (i in 1:n){
      bt <- getbasismatrix(t[i,!is.na(y[i,])&!is.na(t[i,])], basisobj)
      tmp_mtx <- mtx[i,c(1:D)[!is.na(y[i,])&!is.na(t[i,])] + (j-1)*D]
      tmp_mtxx <- rbind(tmp_mtxx, tmp_mtx*bt)
    }
    mtxx <- cbind(mtxx, tmp_mtxx)
  }
}

vy <- vy[!is.na(ty)]
