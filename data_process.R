D <- max(table(dat$id))
unique_id <- unique(dat$id)
tmpt <- tmpy <- matrix(NA,nrow=length(unique_id),ncol=D) 
tmpdatx <- matrix(NA,nrow=length(unique_id),ncol=ncol(datx)) 
tmpmtx <- matrix(NA,nrow=length(unique_id),ncol=ncol(mtx)*D)

for (i in 1:length(unique_id)){
  tmpl <- sum(dat$id==unique_id[i])
  tmpy[i,1:tmpl] <- y[dat$id==unique_id[i]]
  if (tmpl>=2){
    tmpdatx[i,] <- datx[dat$id==unique_id[i],][1,]
  }else{
    tmpdatx[i,] <- datx[dat$id==unique_id[i],]
  }
  tmppmtx <- matrix(NA, nrow=D,ncol=ncol(mtx))
  tmppmtx[1:tmpl,] <- mtx[dat$id==unique_id[i],]
  tmpmtx[i,] <- as.numeric(tmppmtx)
  tmpt[i, 1:tmpl] <- t[dat$id==unique_id[i]]
}
t <- tmpt
y <- tmpy
datx <- tmpdatx
mtx <- tmpmtx