# Add intercept to the design matrix of omics data
newg <- cbind(1,g)

# Define variables of optimization problem corresponding to H_1 under CVXR package
gamma <- Variable(ncol(g)+1,ncls-1)
u <- Variable(n)

# Define the optimization problem corresponding to H_1
tmpw <- w[,1]*DD
obj <- -sum(tmpw*u) 
if (ncls >=2){
  for (cls in 2:ncls){
    tmpw <- w[,cls]*DD
    obj <- obj + sum(tmpw*(cbind(1,g)%*%gamma[,cls-1] - u))
  }
}

# Lasso
if (slambda[1]!=0){
  pen1 <- - slambda[1]*sum(abs(gamma))
  obj <- obj + pen1*n
}

# Group Lasso
if (slambda[2]!=0&length(grpid)==ncol(g)){
  all_grp <- unique(grpid)[unique(grpid)>=0]
  pen2 <- -slambda[2]*norm2(gamma[grpid==all_grp[1],1])
  ngrp <- length(all_grp)
  for (kk in 2:ngrp){
    pen2 <- pen2 -slambda[2]*norm2(gamma[grpid==all_grp[kk],1])
  }
  if (ncls >=3){
    for (pcls in 3:ncls){
      pen2 <- pen2 -slambda[2]*norm2(gamma[grpid==all_grp[1],pcls-1])
      for (kk in 2:ngrp){
        pen2 <- pen2 -slambda[2]*norm2(gamma[grpid==all_grp[kk],pcls-1])
      }
    }
  }
  obj <- obj + pen2*n
}

constraints = list(sum(exp(matrix(newg[1,],nrow=1)%*%gamma-u[1]))+exp(-u[1])<=1)
for (i in 2:n){
  constraints[[i]] = sum(exp(matrix(newg[i,],nrow=1)%*%gamma-u[i]))+exp(-u[i])<=1
}
# Define the optimization under CVXR
prob = Problem(Maximize(obj),constraints)

# Solve the optimization using MOSEK solver
fit1 = psolve(prob,solver='MOSEK')

# Derive the likelihood corresponding to the omics-based sub-model
tmpw <- w[,1]*DD
uhat <- fit1$getValue(u)
ghat <- fit1$getValue(gamma)
H1_loglik <- -sum(tmpw*uhat) 
if (ncls >=2){
  for (cls in 2:ncls){
    tmpw <- w[,cls]*DD
    H1_loglik <- H1_loglik + sum(tmpw*(newg%*%ghat[,cls-1] - uhat))
  }
}
