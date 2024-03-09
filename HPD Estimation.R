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


search <- function(Weights, a)
{
  i = 1
  j = length(Weights)
  mid = 0
  while(j>i+1)
  {
    mid = as.integer((i+j)/2)
    if(Weights[mid]>a){j=mid}
    else{i = mid}
  }
  return(mid)
}

HPD <- function(alpha, Weights)
{
  m = length(Weights)
  vec <- c(NA,NA)
  curr = 10000000
  j1 = (seq(0,as.integer(alpha*m),by = 1))/m
  j2 = j1 + (1-alpha)
  for(i in 1:length(j1))
  {
    theta1 <- Samples[search(Weights, j1[i])]
    theta2 <- Samples[search(Weights, j2[i])]
    if(curr> (theta2-theta1))
    {
      vec = c(theta1,theta2)
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
mu = 7.5
sd = 1.8
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


HPD(0.05,Weights)
qnorm(c(0.25,0.975),mean = posterior_mean,sd = sqrt(posterior_var))

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


HPD(0.05,Weights)
qgamma(c(0.025,0.975),shape = sum(Y)+1, rate = n+1)



