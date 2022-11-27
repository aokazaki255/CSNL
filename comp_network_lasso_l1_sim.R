library(MASS)

#solving constrained Weber problem
constraind_Weber_ADMM <- function(W_hat,R_vec)
{
  n <- nrow(W_hat)
  p <- ncol(W_hat)
  #parameter
  M <- matrix(0,n,p)
  z <- rep(0,p)
  #Lagrangian multipliers
  U <- matrix(0,n,p)
  v <- 0
  #tuning parameters
  rho <- 1
  phi <- 1
  
  Ip<- diag(p)
  Iv<- rep(1,p)
  Iv2<-Iv%*%t(Iv)
  d <- R_vec/rho
  X <- W_hat
  Y <- matrix(0,n,p)
  solved <- solve(rho*n*Ip+phi*Iv2)
  eps <- 100000
  
  while (eps > 0.001)
  {
    r <- rep(0,p)
    old_z <- z 
    for (i in 1:n) 
    {
      Y[i,] <- z - U[i,]/rho
      M[i,] <- prox(d[i],Y[i,],X[i,])
      r <- r + M[i,] + U[i,]/rho
    }
    # z <- r/n
    z <- solved %*% (rho*r-v*Iv)
    for (i in 1:n) 
    {
      U[i,] <- U[i,] + rho*(M[i,]-z)
    }
    v = v + phi*t(Iv)%*%z
    eps <- max(abs(old_z-z))
  }
  return(z)
}

prox <- function(d,y,x)
{
  out  <- y - min(d,sqrt(vec_norm(y-x)))*((y-x)/sqrt(vec_norm(y-x)))
  return(out)
}

#soft_threshold
soft.th = function(x,lambda)
{
  return(sign(x)*pmax(abs(x)-lambda,0))
}

#sum_column-wise for matrix
vec_sum <- function(X)
{
  sumX <- apply(X,2,sum)
  return(sumX)
}

#l2_norm
vec_norm <- function(x)
{
  return(t(x)%*%x)
}

#generating_artificial_data
make_sample = function(n,p,rho_s,sigma)
{
  #mean of multivariate normal distribution
  Mu<- c(rep(log(0.5*p),5),rep(0,p-5))
  
  #variance covariance matrix of multivariate normal distribution
  Sigma <-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      Sigma[i,j] <- rho_s^(abs(i-j))
    }
  }
  
  #generating data W
  W <- mvrnorm(n,Mu,Sigma)
  
  # convert W to compositional data X
  X <- matrix(0,n,p)
  
  for(i in 1:n)
  {
    for(j in 1:p)
    {
      X[i,j] <- exp(W[i,j])/sum(exp(W[i,]))
    }
  }
  
  Z <- log(X)　#taking logarithm for X
  
  #generating epsilon
  Epsilon <- mvrnorm(0,sigma,n=n)
  
  #true regression coefficient vectors
  Beta_ast1 <<- c(1,-0.8,0.6, 0  ,0,-1.5,-0.5,1.2,rep(0,(p-8)))
  Beta_ast2 <<- c(0,-0.5,1  ,1.2,0.1,-1,  0  ,-0.8,rep(0,(p-8)))
  Beta_ast3 <<- c(0,0  ,0,0.8,1,0 ,-0.8,-1,rep(0,(p-8)))
  
  n3  <- n / 3
  
  fir <- n3
  sec <- n3*2
  thi <- n3*3
  
  # calculating y by three clusters
  y1 <<- Z[c(1:fir),]%*%Beta_ast1 + Epsilon[c(1:fir)]
  y2 <<- Z[c((fir+1):sec),]%*%Beta_ast2 + Epsilon[c((fir+1):sec)]
  y3 <<- Z[c((sec+1):thi),]%*%Beta_ast3 + Epsilon[c((sec+1):thi)]
  
  y <- rbind(y1,y2,y3)
  return(list(y,Z)) #return data
}

#making true link information
make_true_R<-function(n)
{
  R <- matrix(0,n,n)
  
  n3  <- n / 3
  
  fir <- n3
  sec <- n3*2
  thi <- n3*3
  
  for (i in 1:fir)
  {
    for(j in 1:fir)
    {
      R[i,j] <- 1 
    }
  }
  
  for (i in (fir+1):sec)
  {
    for(j in (fir+1):sec)
    {
      R[i,j] <- 1 
    }
  }
  
  for (i in (sec+1):thi)
  {
    for(j in (sec+1):thi)
    {
      R[i,j] <- 1 
    }
  }
  
  for (i in 1:n)
  {
    R[i,i] <- 0
  }
  return(R)
}

#true each element of R is obtained in probability P_{R} (=true_p)
make_noised_R <- function(n,true_p)
{
  R <- matrix(0,n,n)
  
  n3  <- n / 3
  
  fir <- n3
  sec <- n3*2
  thi <- n3*3
  c <- 0
  for(i in 1:(fir-1))
  {
    c <- c+1
    r <- fir-c
    for(j in 1:r)
    {
      l <- j+c
      R[i,l] <- sample(x=c(1,0),
                       size=1, 
                       prob=c(true_p,1-true_p)) 
      R[l,i] <- R[i,l]
    }
  }
 
 c <- 0 
  for (i in (fir+1):(sec-1))
  {
    c <- c+1
    r <- sec-c
    for(j in (fir+1):r)
    {
      l <- j+c
      R[i,l] <- sample(x=c(1,0),size=1,prob=c(true_p,1-true_p))
      R[l,i] <- R[i,l]
    }
  }
 
  c <- 0
  for (i in (sec+1):(thi-1))
  {
    c <- c+1
    r <- thi-c
    for(j in (sec+1):r)
    {
      l <- j+c
      R[i,l] <- sample(x=c(1,0),size=1,prob=c(true_p,1-true_p))
      R[l,i] <- R[i,l]
    }
  }
  
  for (i in 1:n)
  {
    R[i,i] <- 0
  }
  return(R)
}

n <- 120
p <- 30
rho_s <- 0.2
sigma <- 0.1
#計画行列の作成
data <- make_sample(n,p,rho_s,sigma)
y <- data[[1]]
Z <- data[[2]]

#link information matrix
true_p = 0.99 #=P_{R} in paper
#R <- make_true_R(n)
R<- make_noised_R(n,true_p)



#matrix consisting of regresion coefficient vectors
W<-matrix(0,n,p)
# W_old <- matrix(0,n,p)
A<-array(rep(0,n),dim=c(n,n,p))　#A[m,l,] =a_{m,l} 
B<-matrix(0,n,p)

#regularization parameters
lambda1 <- 0.5
lambda2 <- 0.2

#Lagrangian multipliers
S<-array(rep(0,n),dim=c(n,n,p))
T<-matrix(0,n,p)
u<-rep(0,n)
#tuning parameters
rho<-1
phi<-1
psi<-1

#pre-caliculation
Ip = diag(p)
Iv = c(rep(1,p))
Iv2 = Iv %*% t(Iv)
#pre-inverse
Wsolve<- vector("list",length = n)
for (i in 1:n)
{
  Wsolve[[i]]  <- solve(2*Z[i,]%*%t(Z[i,])+(rho*(n-1)+phi)*Ip+psi*Iv2)
}


#main-loop
count <- 0
thre <-0.01
loss <- 10000
eps <-100000
while(eps>thre)
{
  count <- count+1
  loss_old <- loss #convergence is judged by value of loss-function
  # W_old  <- W
  #updates of w
  for(i in 1:n)
  {
    W[i,] <- Wsolve[[i]]%*%(2*y[i]*Z[i,]+rho*(vec_sum(A[i,-i,])-vec_sum(S[i,-i,]))-T[i,]+phi*B[i,]-u[i]*Iv)
  }
  #updates of a
  c <- 0
  for(i in 1:(n-1))
  {
    c <- c+1
    r <- n-c
    for(j in 1:r)
    {
      l <- j+c
      theta <- max(1-lambda1*R[i,l]/(rho*sqrt(vec_norm((W[i,]+S[i,l,])-(W[l,]+S[l,i,])))) ,0.5)
      A[i,l,] <- theta*(W[i,]+S[i,l,])+(1-theta)*(W[l,]+S[l,i,])
      A[l,i,] <- (1-theta)*(W[i,]+S[i,l,])+theta*(W[l,]+S[l,i,])
    }
  }
  
  #updates of b
  for (i in 1:n)
  {
    B[i,] <- soft.th((W[i,]+T[i,]/phi),lambda2/phi)
  }
  
  #updates of lagrangian multipliers
  for (i in 1:n)
  {
    for (j in 1:n) 
    {
      if(i!=j)
      {
        S[i,j,] <- S[i,j,]+rho*(W[i,]-A[i,j,])
      }
    }
  }
  T <- T + phi*(W-B)
  u <- u + psi*W%*%Iv
  
  #calculating value ofloss function
  loss <- 0
  for(i in 1:n)
  {
    loss <- loss + (y[i] - t(Z[i,])%*%W[i,])^2
  }
  eps <- abs(loss-loss_old)
  #eps <- max(abs(W-W_old))
}

# return(B)



y_hat <- rep(0,n)

for(i in 1:n)
{
  y_hat[i] <- t(W[i,])%*%Z[i,]
}

res<-y-y_hat

plot(y_hat,res)

# assuming that 1-th data is non observed new i^{ast} data point and calculating corresponding regression coefficient vector w_{i^{ast}} by constrained Weber problem.
W_hat <- B
i_ast <- 1
w_i_ast <- constraind_Weber_ADMM(W_hat[-i_ast,],R[i_ast,-i_ast])
plot(W_hat[i_ast,],w_i_ast)
