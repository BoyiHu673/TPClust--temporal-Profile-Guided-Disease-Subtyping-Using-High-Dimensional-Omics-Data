library(gtools)
library(tidyr)
library(glmnet)
library(CVXR)
library(Rmosek)
library(fda)
library(mvtnorm)
setwd('~/Documents/Columbia_Project_1/submission/')
# Load toy example
dat <- readRDS('toy_example_clinical.RDS')
omics <- readRDS('toy_example_omics.RDS')
g <- omics[[1]]
grpid <- omics[[2]]
grpid <- c(NA,grpid)
# Define the longitudinal response
y <- dat$cognition
# Define the age at each visit
t <- dat$age
# Input the id of individuals
id <- dat$id
unique_id <- unique(id)
# Input the cross-sectional covariates: sex and apoe4
datx <- cbind(dat$sex, dat$apoe4)
# Input the longitudinal covariates: hypertension and diabetes
mtx <- cbind(dat$hypertension, dat$diabetes)
source('data_process.R')

datx <- cbind(1,datx)
age_range <- range(t-min(t,na.rm=T),na.rm=T)[2]
age_min <- min(t,na.rm=T)
# number of time varying coefficients
L <- ncol(datx)
# standardize time points
t <- (t-age_min)/age_range
D <- ncol(t)
# number of longitudinal demographic variables
J <- ncol(mtx)/D
# standardize gene counts
raw_sd <- apply(g,2,sd)
g <- t((t(g) - apply(g,2,mean))/raw_sd)

# group structure of gene expressions
grpid[is.na(grpid)] <- -1
datx <- as.matrix(datx)

# sample size
n <- nrow(g)

# define B-splines
nk <- 8
knots <- quantile(c(0,1),seq(0,1,length=nk))
norder <- 4
nknots <- length(knots)
nb <- nknots + norder - 2
basisobj <- create.bspline.basis(rangeval = c(0,1), nbasis = nb, norder = norder,breaks = knots)
V = eval.penalty(basisobj,int2Lfd(2))
btt <- getbasismatrix(c(1:100)/100,basisobj)

# data vectorization
source('data_vec_real.R')
DD <- apply(!is.na(t+y),1,sum)

org_mtx <- mtx
org_datx <- datx
org_g <- g
org_t <- t
org_y <- y
org_dd <- DD
org_n <- nrow(g)


mtx <- org_mtx
datx <- org_datx
g <- org_g 
t <- org_t
y <- org_y
DD <- org_dd
n <- nrow(g)
source('data_vec_real.R')
# Define number of clusters
ncls <- 3
set.seed(1)
# Initial value for the EM algorithm
shat <- matrix(runif(ncls*2, 0.5, 1),nrow=ncls,ncol=2)
ghat <- cbind(0,matrix(rnorm((ncls-1)*ncol(g)),nrow=ncol(g),ncol=ncls-1))
ghat <- rbind(0,ghat)
that <- matrix(rnorm(ncls*nb*L),nrow=ncls,ncol=L*nb)
rhat <- matrix(rnorm(ncls*nb*J),nrow=ncls,ncol=J*nb)
source('Initialization.R')

# Tuning parameters corresponding to LASSO and group LASSO penalties
tmp_slam <- c(1e-1, 1e-2)
# Tuning parameters corresponding to roughness penalities
tmp_rlam <- rep(1e-6,5)
# Fit TPClust
source('TPClust_EM.R')
ghat[-1,] <- ghat[-1,]/raw_sd
# Identified informative omics
sig_var <- matrix(1,nrow=nrow(ghat),ncol=ncol(ghat)-1)
for (s in 2:ncol(ghat)){
  sig_var[abs(ghat[,s])<=1e-3,s-1] <- 0
}
sig_omics_index <- apply(sig_var[-1,],1,sum)>0
inform_omics <- colnames(g)[sig_omics_index]
print(inform_omics)

# Identified subtypes
subtype1 <- unique_id[w[,1]==apply(w,1,max)]
subtype2 <- unique_id[w[,2]==apply(w,1,max)]
subtype3 <- unique_id[w[,3]==apply(w,1,max)]


# Estimated mean cognitive trajectories across subtypes
ages <- c(1:100)/100*age_range+age_min
plot(ages, btt%*%that[1,1:nb], type='l', ylab='Mean cognition', xlab='Age', main='Mean cognition of Subtype 1')
plot(ages, btt%*%that[2,1:nb], type='l', ylab='Mean cognition', xlab='Age', main='Mean cognition of Subtype 2')
plot(ages, btt%*%that[3,1:nb], type='l', ylab='Mean cognition', xlab='Age', main='Mean cognition of Subtype 3')

# Estimated time-varying effects of sex across subtypes
ages <- c(1:100)/100*age_range+age_min
plot(ages, btt%*%that[1,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of sex for Subtype 1')
plot(ages, btt%*%that[2,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of sex for Subtype 2')
plot(ages, btt%*%that[3,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of sex for Subtype 3')

# Estimated time-varying effects of APOE4 across subtypes
ages <- c(1:100)/100*age_range+age_min
plot(ages, btt%*%that[1,1:nb+nb*2], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 1')
plot(ages, btt%*%that[2,1:nb+nb*2], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 2')
plot(ages, btt%*%that[3,1:nb+nb*2], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 3')


# Estimated time-varying effects of hypertension across subtypes
ages <- c(1:100)/100*age_range+age_min
plot(ages, btt%*%rhat[1,1:nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 1')
plot(ages, btt%*%rhat[2,1:nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 2')
plot(ages, btt%*%rhat[3,1:nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 3')

# Estimated time-varying effects of diabetes across subtypes
ages <- c(1:100)/100*age_range+age_min
plot(ages, btt%*%rhat[1,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 1')
plot(ages, btt%*%rhat[2,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 2')
plot(ages, btt%*%rhat[3,1:nb+nb], type='l', ylab='Effect on cognition', xlab='Age', main='Time-varying effect of APOE4 for Subtype 3')
