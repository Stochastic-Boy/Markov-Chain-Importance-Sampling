set.seed(100)

## Density function
dens1 <- function(x){return(exp(-x^2)*(2+sin(5*x)+sin(2*x)))}

## Gradient of Density function
grad_dens1 <- function(x){return(-2*x+ (5*cos(5*x)+2*cos(2*x))/(2+sin(5*x)+sin(2*x)))}


## Metropolis - Hastings Algorithm for given density function
mh_sampler <- function(dens, start = 1, nreps = 1000, prop_sd = 1)
{
    theta <- rep(0,nreps)
    theta[1] <- start
    prop<- rep(0,nreps)
    prop <- start
    
    for (i in 2:nreps)
    {
      theta_star <- rnorm(1, mean = theta[i-1], sd = prop_sd)
      alpha = dens(theta_star)/dens(theta[i-1])
      U = runif(1)
      prop[i] = theta_star
      if(U<alpha){theta[i] <- theta_star}else{theta[i] <- theta[i-1]}
    }
    return(cbind(theta,prop))
}


## Unadjusted Langevin algorithm for given density function
ula_sampler <- function(dens, start = 2, nreps = 1000, paramet = 0.7)
{
  theta <- rep(0,nreps)
  theta[1] <- start
  for (i in 2:nreps)
  {
    theta_star <- rnorm(1, mean = theta[i-1]+paramet*grad_dens(theta[i-1]), sd = 2*paramet)
    theta[i] <- theta_star
  }
  return(cbind(theta,theta))
}


## Weight function for the given data
weight <- function(y,data){return(dens(y)/rhoymh(y,data))}
rhoymh <- function(y,data)
{
  k = dim(data)[1]
  ans = 0
  ans = sum(dnorm(y,mean=data[,1],sd=rep(1,k)))
  return(ans/k)
}
  

## Estimation of function fun 
estimate <- function(fun,data,weight)
{
  sum1 = 0
  sum2 = 0
  k = dim(data)[1]
  for(i in 1:k)
  {
    sum1 = sum1 + weight(data[i,2],data)*fun(data[i,2])
    sum2 = sum2 + weight(data[i,2],data) 
  }
  return(sum1/sum2)
}


## Some general function
fun1 <- function(x){return(x)}
fun2 <- function(x){return(x^2)}


##########################************** Implementation ***********************#####################################



## Generate and plot data from MHA Sampler
data1 <- mh_sampler(dens = dens, start=1, nreps = 10000)
plot(density(data[,1]))

## Generate and plot data from ULA Sampler
data2 <- ula_sampler(dens = dens, start=1, nreps = 10000)
plot(density(data))

estimate(fun1,data1,weight)
estimate(fun2,data1,weight)

dat2 <- vector(length=dim(data1)[1])
for(i in 1:length(dat1))
{
  dat2[i] = rhoymh(data1[i,2],data1)
}

