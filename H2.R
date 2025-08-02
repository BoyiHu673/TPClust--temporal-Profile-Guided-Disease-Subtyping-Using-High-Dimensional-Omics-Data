# Define the design matrix for H_{2,k}
tmpvw <- sqrt(vw)
wy <- tmpvw*vy
waa <- tmpvw*aa
waavar <- vw*aavar
wxx <- tmpvw*xx

if (!is.na(J)) wxx <- cbind(wxx,tmpvw*mtxx)

# Define the roughness penalties for H_{2,k}
tmpm <- rlambda[1]*V*n
if (L>=2){
  for (jj in 2:L){
    tmpm <- bdiag(tmpm,rlambda[jj]*V*n)
  }
}

if (!is.na(J)){
  tmpm <- bdiag(tmpm,rlambda[L+1]*V*n)
  if (J>=2){
    for (jj in 2:J){
      tmpm <- bdiag(tmpm,rlambda[jj+L]*V*n)
    }
  }
} 

# Find the maximizer of H_{2,k}
h2hat <- solve(t(wxx)%*%wxx + tmpm)%*%t(wxx)%*%as.matrix(wy-waa)
tmp_sse <- sum((wy-waa-wxx%*%h2hat)^2)
tmp_sse <- tmp_sse + sum(w[,k]*DD*alpha_var[,k])

# Derive the estimates for parameters
shat[k,1] <- (as.numeric(tmp_sse)/sum(w[,k]*DD))^0.5
shat[k,2] <- ((sum(w[,k]*DD*alpha_var[,k]+w[,k]*DD*alpha[,k]^2))/sum(w[,k]*DD))^0.5
that[k,] <- h2hat[1:(nb*L)]
if (!is.na(J)) rhat[k,] <- h2hat[1:(nb*J)+nb*L]