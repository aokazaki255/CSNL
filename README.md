# CSNL
Authors: Akira Okazaki and Shuichi Kawano
## Introduction
This is an implementation of proposed method in "Multi-task Learning for Compositional Data via Sparse Network Lasso".

1. Import some libraries and our source codes.
  ```
  library(MASS)
  source("R/make_dataset.R")
  source("R/cal_ADMM.R")
  ```
  
2. Generate an artificial data
  Set the following parameters, and generate a dataset.
    - n: Sample size (choose a multiple of three)
    - p: The dimension of explanatory variables
    - rho_s: Correlation of explanatory variables
    - sigma: Standard deviation of error term
    - true_p: Probability of true observation for R
    
   ```
   n <- 120
   p <- 30
   rho_s <- 0.2
   sigma <- 0.1
   
   data <- make_sample(n,p,rho_s,sigma)
   y <- data[[1]]
   Z <- data[[2]]
      
   true_p = 0.99 #=P_{R} in our paper
   R<- make_noised_R(n,true_p)
   ```
  
3. perform the proposed method

    - estimation of proposed method 
       ```
       W_hat <- CSNL_estimator(y,Z,R,lambda1,lambda2)
       ``` 
    - estimation of constrained Weber problem
       ```
       i_ast <- 1
       w_i_ast <- constraind_Weber_ADMM(W_hat[-i_ast,],R[i_ast,-i_ast])
       ```

### Licence

These code are free, non-commercial and open source.

