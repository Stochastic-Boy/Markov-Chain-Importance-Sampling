library(mvtnorm)

compute_gradient <- function(x) {
  # Gradient of the log-density of the target distribution (Gaussian with mean 0 and covariance matrix identity)
  return(-x)
}

multivariate_ULA_sampler <- function(initial_state, step_size, num_samples, mean) {
  # Function to perform multivariate Unadjusted Langevin Algorithm (ULA) sampling
  
  # Initialize variables
  d <- length(initial_state)
  samples <- matrix(0, nrow=num_samples, ncol=d)
  x <- initial_state
  
  # Perform sampling
  for (i in 1:num_samples)
  {
    gradient <- compute_gradient(x)
    x <- rmvnorm(1, mean = x - epsilon*(x-mean), sigma = 2 * epsilon*diag(d))
    samples[i,] <- x
  }
  
  return(samples)
}

# Example usage:
set.seed(123)
initial_state <- rep(-1.5, 5)  # Initial state
mean = rep(0,5)
step_size <- 0.001  # Step size
num_samples <- 100000  # Number of samples

samples <- multivariate_ULA_sampler(initial_state, step_size, num_samples, mean)
print(head(samples))  # Print first few samples

samp = sort(samp)
g1 <- density(samp)
weights <- vector(length = length(samp))
weights = dnorm(samp, mean = 0, sd = 1)/(approx(g1$x,g1$y,xout=samp)$y)
weights = weights/sum(weights)
Weights = cumsum(weights)

HPD <- function(alpha, Weights, Samples)
{
  m = length(Weights)
  vec <- c(NA,NA)
  curr = 10000000
  j1 = (seq(1,floor((alpha)*m),by = 1))/m
  j2 = j1 + ((1-alpha))
  for(i in 1:(length(j1)-1))
  {
    theta <- c(Samples[which(Weights>j1[i])[1]], Samples[which(Weights >= j2[i])[1]])
    if(curr > (theta[2]-theta[1]))
    {
      vec = theta
      curr = (theta[2]-theta[1])
    }
  }
  return(vec)
}

HPD(0.05,Weights, samp)
qnorm(c(0.025,0.975),mean = 0,sd = 1)


alpha <- c(0.25,0.5,0.75)
MSE <- rep(0,3)
qu <- qgamma(c(0.25,0.5,0.75),shape = sum(Y)+1, rate = n+1)
for(i in 1:10000)
{
  samples <- multivariate_ULA_sampler(initial_state, step_size, num_samples, mean)
  samp <- samples[,1]
  samp = sort(samp)
  g1 <- density(samp)
  weights <- vector(length = length(samp))
  weights = dnorm(samp, mean = 0, sd = 1)/(approx(g1$x,g1$y,xout=samp)$y)
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

