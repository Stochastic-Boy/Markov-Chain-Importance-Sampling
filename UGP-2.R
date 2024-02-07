set.seed(100)

## Density function
dens <- function(x){return(exp(-x^2)*(2+sin(5*x)+sin(2*x)))}

## Gradient of Density function
grad_dens <- function(x){return(-2*x+ (5*cos(5*x)+2*cos(2*x))/(2+sin(5*x)+sin(2*x)))}


## Metropolis - Hastings Algorithm for given density function
mh_sampler <- function(dens, start = 1, nreps = 1000, prop_sd = 1)
{
  theta <- rep(0,nreps)
  theta[1] <- start
  prop<- rep(0,nreps)
  
  for (i in 2:nreps)
  {
    theta_star <- rnorm(1, mean = theta[i-1], sd = prop_sd)
    alpha = dens(theta_star)/dens(theta[i-1])
    U = runif(1)
    prop[i-1] = theta_star
    if(U<alpha){theta[i] <- theta_star}else{theta[i] <- theta[i-1]}
  }
  prop[nreps] <- rnorm(1, mean = theta[nreps], sd = prop_sd)
  return(cbind(theta,prop))
}


## Unadjusted Langevin algorithm for given density function
ula_sampler <- function(dens, start = 2, nreps = 1000, paramet = 0.7)
{
  theta <- rep(0,nreps)
  theta[1] <- start
  prop <- rep(0,nreps)
  # q_y <- rep(0,nreps)
  for (i in 2:nreps)
  {
    theta_star <- rnorm(1, mean = theta[i-1]+paramet*grad_dens(theta[i-1]), sd = sqrt(2*paramet))
    #  q_y[i-1] <- dnorm(theta_star, mean = theta[i-1]+paramet*grad_dens(theta[i-1]), sd = sqrt(2*paramet))
    prop[i-1] <- theta_star
    theta[i] <- theta_star
  }
  prop[nreps] <- rnorm(1, mean = theta[nreps]+paramet*grad_dens(theta[nreps]), sd = sqrt(2*paramet))
  # q_y[nreps] <- dnorm(prop[nreps], mean = theta[nreps]+paramet*grad_dens(theta[nreps]), sd = sqrt(2*paramet))
  return(cbind(theta,prop))
}



rhoymh <- function(y,data)
{
  k = dim(data)[1]
  ans = 0
  ans = sum(dnorm(y,mean=data[,1],sd=rep(1,k)))
  return(ans/k)
}

rhoyula <- function(y,data, paramet)
{
  k = dim(data)[1]
  ans = 0
  ans = mean(dnorm(y,mean = data[,1]+paramet*grad_dens(data[,1]),sd= sqrt(2*paramet)))
  return(ans)
}


## Weight function for the given data
weight_mh <- function(y,data, paramet = .01){return(dens(y)/rhoymh(y,data))}
weight_ula <- function(y,data, paramet){return(dens(y)/ rhoyula(y, data, paramet) )}

## Estimation of function fun 
estimate <- function(fun,data,weight, paramet = .01)
{
  sum1 = 0
  sum2 = 0
  k = dim(data)[1]
  for(i in 1:k)
  {
    sum1 = sum1 + weight(data[i,2],data, paramet)*fun(data[i,2])
    sum2 = sum2 + weight(data[i,2],data, paramet) 
  }
  return(sum1/sum2)
}


## Some general function
fun1 <- function(x){return(x)}
fun2 <- function(x){return(x^2)}


##########################************** Implementation ***********************#####################################



## Generate and plot data from MH Sampler

data1 <- mh_sampler(dens = dens, start=1, nreps = 1e4)
plot(density(data1[,1]))
lines(density(data1[,2]), col = "tomato")

y_dum <- seq(-4, 4, length = 100)
dat2 <- numeric(length(y_dum))
for(i in 1:length(dat2))
{
  dat2[i] = rhoymh(y_dum[i],data1)
}
lines(y_dum, dat2, col = "blue")


paramet = .2
## Generate and plot data from ULA Sampler
data2 <- ula_sampler(dens = dens, start=1, nreps = 1e4, paramet = paramet)
lines(density(data2), col = "red")

y_dum <- seq(-4, 4, length = 100)
dat2 <- numeric(length(y_dum))
for(i in 1:length(dat2))
{
  dat2[i] = rhoyula(y_dum[i],data2, paramet = paramet)
}
lines(y_dum, dat2, col = "blue")



estimate(fun1,data1,weight_mh)
mean(data1[, 1])
estimate(fun2,data1,weight_mh)
mean(data1[, 1]^2)

estimate(fun1,data2,weight_ula, paramet = paramet)
mean(data1[, 1])
estimate(fun2,data2,weight_ula, paramet = paramet)
mean(data1[, 1]^2)





