
#main functions

#rho,phi,psi are tuning parameters related to convergence of ADMM
#thre is convergence parameter we recommend a small value as possible
#estimating proposed method
CSNL_estimator <- function(y,Z,R,lambda1,lambda2,rho=1,phi=1,psi=1,thre=0.001)
{
  n=nrow(Z)
  p=ncol(Z)
  
  #matrix consisting of regresion coefficient vectors
  W<-matrix(0,n,p)
  # W_old <- matrix(0,n,p)
  A<-array(rep(0,n),dim=c(n,n,p))　#A[m,l,] =a_{m,l} 
  B<-matrix(0,n,p)
  
  #Lagrangian multipliers
  S<-array(rep(0,n),dim=c(n,n,p))
  T<-matrix(0,n,p)
  u<-rep(0,n)
  
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
  
  return(B)
}

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


##other functions

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
}#
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
