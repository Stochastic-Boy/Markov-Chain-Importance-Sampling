library(mvtnorm)

# Function to generate samples using Unadjusted Langevin Algorithm for Gaussian target
ula_gaussian <- function(n_samples, epsilon, mean, variance, start) 
{
  samples <- numeric(n_samples)
  x <- start
  for (i in 1:n_samples) 
  {
    x = rnorm(1, mean= x - epsilon*(x-mean)/variance, sd = sqrt(2 * epsilon))
    samples[i] <- x
  }
  return(samples)
}

# Function to generate samples using Unadjusted Langevin Algorithm for Gamma target
ula_gamma <- function(n_samples, epsilon, shape, rate) 
{
  samples <- numeric(n_samples)
  x <- rgamma(1, shape = shape, rate= rate)
  for (i in 1:n_samples)
  {
    x = rnorm(1, mean = x + epsilon * ((shape-1)/x - rate), sd = sqrt(2 * epsilon))
    samples[i] <- x
  }
  return(samples)
}



HPD <- function(alpha, Weights, Samples)
{
  m = length(Weights)
  vec <- c(NA,NA)
  curr = 10000000
  j1 = (seq(1,floor((alpha)*m),by = 1))/m
  j2 = j1 + ((1-alpha))
  for(i in 1:length(j1))
  {
    theta <- c(Samples[which(Weights>j1[i])[1]],Samples[which(Weights>=j2[i])[1]])
    if(curr > (theta[2]-theta[1]))
    {
      vec = theta
      curr = (theta[2]-theta[1])
    }
  }
  return(vec)
}


# Example usage for Gaussian target
set.seed(123)
n_samples <- 1000
epsilon <- 0.1
samples_gaussian <- ula_gaussian(n_samples, epsilon,0,1,1)
plot(density(samples_gaussian))


# Example usage for Gamma target
set.seed(1)
n_samples <- 1000
epsilon <- 0.01
shape <- 2
rate <- 1
samples_gamma <- ula_gamma(n_samples, epsilon, shape, rate)
plot(density(samples_gamma))
 

#######################################################################################
#######################################################################################

# Let, in the case of gaussian ULA sampler, we need to find the quantiles for mean of the 
# underlying distribution when population variance known.

# prior mu ~ N(0,100)
set.seed(4452)
mu = 0
sd = 1
n = 50
Y <- rnorm(n, mean = mu, sd = sd)
posterior_mean <- (sum(Y)/sd^2)/((n/sd^2) + (1/100))
posterior_var <- 1/((n/sd^2) + (1/100))

m = 100000
Samples <- ula_gaussian(m,epsilon = 0.0001,posterior_mean,posterior_var,start = 0.000001)
plot(density(Samples))
Samples <- sort(Samples)
g <- density(Samples)
weights <- vector(length = m)
weights = dnorm(Samples,mean = posterior_mean,sd = sqrt(posterior_var))/(approx(g$x,g$y,xout=Samples)$y)
weights = weights/sum(weights)
Weights = cumsum(weights)


HPD(0.05,Weights, Samples)
qnorm(c(0.025,0.975),mean = posterior_mean,sd = sqrt(posterior_var))

alpha <- c(0.25,0.5,0.75)
MSE <- rep(0,3)
qu <- qnorm(c(0.25,0.5,0.75),mean = posterior_mean,sd = sqrt(posterior_var))
for(i in 1:1000)
{
  Samples <- ula_gaussian(m,epsilon = 0.0001,posterior_mean,posterior_var,start = 0.000001)
  Samples <- sort(Samples)
  g <- density(Samples)
  weights <- vector(length = m)
  weights = dnorm(Samples,mean = posterior_mean,sd = sqrt(posterior_var))/(approx(g$x,g$y,xout=Samples)$y)
  weights = weights/sum(weights)
  Weights = cumsum(weights)
  q1 <- Samples[which(Weights>alpha[1])]
  q2 <- Samples[which(Weights>alpha[2])]
  q3 <- Samples[which(Weights>alpha[3])]
  MSE[1] <- MSE[1] + (q1-qu[1])^2
  MSE[2] <- MSE[2] + (q2-qu[2])^2
  MSE[3] <- MSE[2] + (q3-qu[3])^2
}

MSE <- MSE/1000
MSE
## >>> 2.361802e-07 3.490802e-08 3.491659e-08

MSE <- rep(0,3)
m = 100000
qu <- qnorm(c(0.25,0.5,0.75),mean = posterior_mean,sd = sqrt(posterior_var))
for(i in 1:10000)
{
  Samples <- ula_gaussian(m,epsilon = 0.0001,posterior_mean,posterior_var,start = 0.000001)
  Samples <- sort(Samples)
  q1 <- Samples[alpha[1]*m]
  q2 <- Samples[alpha[2]*m]
  q3 <- Samples[alpha[3]*m]
  MSE[1] <- MSE[1] + (q1-qu[1])^2
  MSE[2] <- MSE[2] + (q2-qu[2])^2
  MSE[3] <- MSE[3] + (q3-qu[3])^2
}

MSE <- MSE/10000
MSE
## >>> 9.630869e-05 8.786842e-05 9.504091e-05
#########################################################################################
#########################################################################################


## Now Let try this approach for poison distribution.

lambda = 14
n = 50
Y = rpois(n,lambda)

# Prior: Gamma(1,1)
# Posterior: Gamma(sum(Y)+1, n+1)

m = 100000
Samples <- ula_gamma(m,epsilon = 0.001,shape = sum(Y)+1, rate = (n+1))
plot(density(Samples))
Samples <- sort(Samples)
g <- density(Samples)
weights <- vector(length = m)
weights = dgamma(Samples,shape = sum(Y)+1, rate = n+1)/(approx(g$x,g$y,xout=Samples)$y)
weights = weights/sum(weights)
Weights = cumsum(weights)

fun<-function(x)
{
  k <- qgamma(c(x+0.025,x+0.975),shape = sum(Y)+1, rate = n+1)
  return(k[2]-k[1])
}

x <- optimize(fun,lower=0,upper = 0.025)$minimum
HPD_exact <- qgamma(c(x+0.025,x+0.975),shape = sum(Y)+1, rate = n+1)
HPD(0.05,Weights,Samples)


alpha <- c(0.25,0.5,0.75)
MSE <- rep(0,3)
qu <- qgamma(c(0.25,0.5,0.75),shape = sum(Y)+1, rate = n+1)
for(i in 1:10000)
{
  Samples <- ula_gamma(m,epsilon = 0.001,shape = sum(Y)+1, rate = (n+1))
  Samples <- sort(Samples)
  g <- density(Samples)
  weights <- vector(length = m)
  weights = dgamma(Samples,shape = sum(Y)+1, rate = n+1)/(approx(g$x,g$y,xout=Samples)$y)
  weights = weights/sum(weights)
  Weights = cumsum(weights)
  q1 <- Samples[which(Weights>alpha[1])]
  q2 <- Samples[which(Weights>alpha[2])]
  q3 <- Samples[which(Weights>alpha[3])]
  MSE[1] <- MSE[1] + (q1-qu[1])^2
  MSE[2] <- MSE[2] + (q2-qu[2])^2
  MSE[3] <- MSE[3] + (q3-qu[3])^2
}

MSE <- MSE/10000
MSE

## >>> 3.686399e-06 6.695139e-07 4.705776e-06


MSE <- rep(0,3)
m = 100000
qu <- qgamma(c(0.25,0.5,0.75),shape = sum(Y)+1, rate = n+1)
for(i in 1:10000)
{
  Samples <- ula_gamma(m, epsilon = 0.001,shape = sum(Y)+1, rate = (n+1))
  Samples <- sort(Samples)
  q1 <- Samples[alpha[1]*m]
  q2 <- Samples[alpha[2]*m]
  q3 <- Samples[alpha[3]*m]
  MSE[1] <- MSE[1] + (q1-qu[1])^2
  MSE[2] <- MSE[2] + (q2-qu[2])^2
  MSE[3] <- MSE[3] + (q3-qu[3])^2
}

MSE <- MSE/10000
MSE
## >>> 0.001605519 0.001539897 0.001783654