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
- **TPClust_EM.R**: This is the main script that implements the EM algorithm for fitting the TPClust framework.
- **data_processing.R**: This script handles the data preprocessing and transforms the input data into the proper format for TPClust.
- **data_vec.R**: The script is used for vectorizing the data, preparing it in a form that can be input into the TPClust algorithm.
- **H1_sparse.R**: This script fits a weighted multinomial logistic regression model with Lasso and group Lasso penalties, specifically applied to omics data within the EM algorithm, where the weights are derived from the E-step of the EM algorithm.
- **H2.R**: This script fits the weighted outcome association model with roughness penalties within the EM algorithm, where the weights are derived from the E-step of the EM algorithm.
- **Initialization.R**: This script generates the initial points for the EM algorithm.
- **example.R**: This is a demonstration script that runs the TPClust framework using a toy example. It's useful for testing the functionality of the code and understanding how TPClust works in practice.
