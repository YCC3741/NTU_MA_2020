rm(list=ls())                      # clear memory

# Metropolis Hastings

set.seed(1);x <- rnorm(20)*2+1 ## simulated data 
n.rep <- 10000; n.accept <- 0

theta <- matrix(0,2,n.rep) ## storage for sim. values 
ll0 <- sum(dnorm(x,mean=theta[1,1],sd=exp(theta[2,1]),log=TRUE)) 

for (i in 2:n.rep) { ## The MH loop
  theta[,i] <- theta[,i-1] + rt(2,df=3)*.5 ## proposal 
  ll1 <- sum(dnorm(x,mean=theta[1,i],
  sd=exp(theta[2,i]),log=TRUE))
  if (exp(ll1-ll0)>runif(1)) { ## MH accept/reject
  ll0 <- ll1; n.accept <- n.accept + 1 ## accept 
  } else theta[,i] <- theta[,i-1] ## reject
}

n.accept/n.rep ## proportion of proposals accepted


layout(matrix(c(1,2,1,2,3,4),2,3))
plot(1:n.rep,theta[1,],type="l",xlab="iteration",
     ylab=expression(mu))
plot(1:n.rep,exp(theta[2,]),type="l",xlab="iteration",
     ylab=expression(sigma))
hist(theta[1,-(1:1000)],main="",xlab=expression(mu))
hist(exp(theta[2,-(1:1000)]),main="",
     xlab=expression(sigma))


# Gibbs sampling

n <- 20;set.seed(1);x <- rnorm(n)*2+1 ## simulated data 
n.rep <- 10000;
thetag <- matrix(0,2,n.rep)
a <- 1; b <- .1; c <- 0; d <- 100 ## prior constants 
xbar <- mean(x) ## store mean
thetag[,1] <- c(mu <- 0,phi <- 1) ## initial guesses 

for (j in 2:n.rep) { ## the Gibbs sampling loop
  mu <- rnorm(1,mean=(d*n*xbar+phi*c)/(d*n+phi),
            sd=sqrt(phi*d/(d*n+phi)))
  phi <- 1/rgamma(1,n/2+a,sum((x-mu)^2)/2+b)
  thetag[,j] <- c(mu,phi) ## store results 
}


layout(matrix(c(1,2,1,2,3,4),2,3))
plot(1:n.rep,thetag[1,],type="l",xlab="iteration",
     ylab=expression(mu))
plot(1:n.rep,thetag[2,],type="l",xlab="iteration",
     ylab=expression(sigma))
hist(thetag[1,-(1:1000)],main="",xlab=expression(mu))
hist(thetag[2,-(1:1000)],main="",
     xlab=expression(sigma))


# Checking for convergence

qtplot <- function(theta,n.plot=100,ylab="") { ## simple MCMC chain diagnostic plot 
  cuq <- Vectorize(function(n,x) ## cumul. quantile func.
    as.numeric(quantile(x[1:n],c(.025,.5,.975))),
    vectorize.args="n")
  n.rep <- length(theta)
  plot(1:n.rep,theta,col="lightgrey",xlab="iter", ylab=ylab,type="l")
  iter <- round(seq(1,n.rep,length=n.plot+1)[-1])
  tq <- cuq(iter,theta)
  lines(iter,tq[2,])
  lines(iter,tq[1,],lty=2);lines(iter,tq[3,],lty=2)
}

layout(matrix(c(1,2,1,2,3,4),2,3))
qtplot(theta[1,],ylab=expression(mu)) 
qtplot(exp(theta[2,]),ylab=expression(sigma)) 
acf(theta[1,])
acf(exp(theta[2,]))
