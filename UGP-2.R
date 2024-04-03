set.seed(100)


dens <- function(x){return(exp(-x^2)*(2+sin(5*x)+sin(2*x)))}
grad_dens <- function(x){return(-2*x+ (5*cos(5*x)+2*cos(2*x))/(2+sin(5*x)+sin(2*x)))}

## A bi-variate Density Function
ring2D <- function(x){ exp(-5*(x[1]^2+x[2]^2-1)^2)}
grad_ring2D <- function(x){return(-20*(x[1]^2+x[2]^2-1)*(x[1]+x[2]))}

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
  for (i in 2:nreps)
  {
    theta_star <- rnorm(1, mean = theta[i-1]+paramet*grad_dens(theta[i-1]), sd = sqrt(2*paramet))
    prop[i-1] <- theta_star
    theta[i] <- theta_star
  }
  prop[nreps] <- rnorm(1, mean = theta[nreps]+paramet*grad_dens(theta[nreps]), sd = sqrt(2*paramet))
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

estimate2 <- function(fun,data, paramet = .01)
{
  sum1 = 0
  sum2 = 0
  k = dim(data)[1]
  g = density(data[,2])
  y = data[,2]
  for(i in 1:k)
  {
    sum1 = sum1 + (dens(y[i])/approx(g$x,g$y,xout = y[i])$y)*fun(data[i,2])
    sum2 = sum2 + dens(y[i])/approx(g$x,g$y,xout = y[i])$y
  }
  return(sum1/sum2)
}


## Some general univariate function
fun1 <- function(x){return(x)}
fun2 <- function(x){return(x^2)}

## Multivariate Function

testfunction1 <- function(x){d = length(x); return(sum(x^3)/d)}

testfunction2 <- function(x){d = length(x); return(sum(x^2)/d)}

testfunction2 <- function(x){d = length(x); return(sum(x)/d)}

testfunction4 <- function(x){d = length(x); return(sum(exp(x))/d)}


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


paramet = 0.2
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
estimate2(fun1,data1)
mean(data1[, 1])
estimate(fun2,data1,weight_mh)
mean(data1[, 1]^2)

estimate(fun1,data2,weight_ula, paramet = paramet)
mean(data1[, 1])
estimate(fun2,data2,weight_ula, paramet = paramet)
mean(data1[, 1]^2)
################################# Variance of Estimators ########################################

est_mh <- vector(length=1000)

for(i in 1:1000)
{
  data_mh <- mh_sampler(dens = dens, start=1, nreps = 1e3)
  est_mh[i] <- estimate(fun1,data_mh,weight_mh)
}

var(est_mh)
# 0.003161821

est_mh_kde <- vector(length=1000)
for(i in 1:1000)
{
  data_mh <- mh_sampler(dens = dens, start=1, nreps = 1e4)
  est_mh_kde[i] <- estimate2(fun1,data_mh)
}
var(est_mh_kde)
#   7.337105e-07

mean(est_mh)
# 0.1868874
mean(est_mh_kde)
# 0.1880719

data_mh <- mh_sampler(dens = dens, start=1, nreps = 1e6)
mean(data_mh)
# 0.1899867


data<- matrix(NA, ncol=2,nrow=7)
j = 1
for(i in c(100,500,1000,5000,10000,50000,100000))
{
  data_mh <- mh_sampler(dens = dens, start=1, nreps = i)
  data[j,1] = system.time({out = estimate(fun1,data_mh,weight_mh)})[1]
  data[j,2] = system.time({out = estimate2(fun1,data_mh)})[1]
  j = j +1
}

x = c(100,500,1000,5000,10000,50000,100000)
par(mfrow = c(1,2))
plot(x = x, y = data[,1], "l",col = "red",lwd="2",main="Performance of MCIS", xlab = "Number of Samples",ylab = "System time")
plot(x = x, y = data[,2], "l",col = "blue",lwd="2",main="Performace of KDE estimate", xlab = "Number of Samples",ylab = "System time")

