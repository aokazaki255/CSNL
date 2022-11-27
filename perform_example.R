library(MASS)
source("R/make_dataset.R")
source("R/cal_ADMM.R")


#settings for simulation
# n: sample-size
# p: dimensionality
# rho_s: correlation between features 
# sigma: standard deviation of error term
n <- 120
p <- 30
rho_s <- 0.2
sigma <- 0.1


#generating artificial data
data <- make_sample(n,p,rho_s,sigma)
y <- data[[1]]
Z <- data[[2]]


#link information matrix

#R <- make_true_R(n)

true_p = 0.99 #=P_{R} in paper
R<- make_noised_R(n,true_p) 


#regularization parameters
lambda1 <- 2
lambda2 <- 0.4

#estimation of proposed method
W_hat <- CSNL_estimator(y,Z,R,lambda1,lambda2)


##example for calculating w_{i}^{ast}
# assuming that 1-th data is non observed new i^{ast} data point and calculating corresponding regression coefficient vector w_{i^{ast}} by constrained Weber problem.
i_ast <- 1
w_i_ast <- constraind_Weber_ADMM(W_hat[-i_ast,],R[i_ast,-i_ast])
