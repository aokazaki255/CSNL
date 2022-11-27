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
  Epsilon <- rnorm(0,sigma,n=n)
  
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
