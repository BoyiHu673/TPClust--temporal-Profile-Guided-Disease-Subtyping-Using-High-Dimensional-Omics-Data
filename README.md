# TPClust--temporal-Profile-Guided-Disease-Subtyping-Using-High-Dimensional-Omics-Data
## This markdown file contains R code for the manuscript entitled "TPClust: Temporal Profile-Guided Subtyping Using High-Dimensional Omics Data"

## Our implementation uses the MOSEK solver for optimization problems with LASSO and group LASSO penalties. Therefore, the installation of MOSEK solver in R is required.
---------------------------------------------------------------
1. Download the correct version of MOSEK based on your machine from the official website (https://www.mosek.com/downloads/)
2. Set up a license following the instructions (https://docs.mosek.com/11.0/licensing/index.html). For academic faculty, students or staff, Academic Licenses for MOSEK can be request from https://www.mosek.com/products/academic-licenses/
3. Follow installation instructions (https://docs.mosek.com/11.0/rmosek/install-interface.html) to install MOSEK in R.

---------------------------------------------------------------
## This repository contains the R code used to excute TPClust method and one illustration of fitting TPClust method on a toy example. The descriptions of R files in this repository are given below
---------------------------------------------------------------
- TPClust_EM: The R code that perfomrs the EM algorithm used to fit TPClust framework
- data_processing.R: The R code used to transfer the input data into the right format for TPClust
- data_vec.R: The R code used to vectorize and prepare the data
- H1_sparse.R: Fit the weighted multinomial logistic regression with Lasso and group Lasso penaltie corresponding to omics data within the EM algorithm
- H2.R: Fit the weighted outcome association model with roughnes penalties witin the EM algorithm
- Initialization.R: The R code used to generate the initial point of the EM algorithm
- example.R: The R code excute a demonstration for TPClust_EM using a toy example
