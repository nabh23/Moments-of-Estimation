library(maxLik)


set.seed(1001)

N <- 100

x <- rnorm(N, mean = 3, sd = 2)

mean(x)   # 2.9983
sd(x)     # 2.2889


# log likelihood function for normal distribution

LL <- function(mu, sigma) {
  R <- dnorm(x, mu, sigma)
 # print(R)
  -sum(log(R))
  #print(R)
  #print (sigma) 
}


mle(LL, start = list(mu = 1, sigma =1))

# checking for negatives in standard deviation
dnorm(x, 1, -1)

# to fix this issue, applying constraints on the parameters   # L-BFGS-B alows box constraints
mle(LL, start = list(mu = 1, sigma=1), method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(Inf,Inf))



# changing the start values
mle(LL, start = list(mu = 0, sigma =0))


#######################

LL_maxLik <- function(mu, sigma, data) {
  R <- dnorm(x, mu, sigma)
  # print(R)
  logLik = -sum(log(R))
  #print(R)
  #print (sigma) 
}


LL_maxLik <- function(mu, sigma, data) {
  R <- dnorm(data, mu, sigma)
  # print(R)
  logLik = -sum(log(R))
  #print(R)
  #print (sigma) 
}

maxLik(LL_maxLik, start = list(mu = 1,sigma = 1), data = data)

maxLik(LL_maxLik, start = c(mu = 1,sigma = 1))






##########################################################
# Estimating 2 parameter models
# Normal Distribiution


set.seed(1001)

N <- 100
# Generating normally distribiuted input vector
x <- rnorm(N, mean = 3, sd = 2)

n = length(x)

# log-likelihood function for Normal distribution

loglik_normal <- function(theta,x){
  # theta : parameters c(mu, sigma2)
  
  mu <- theta[1]
  sig2 <- theta[2]
  n = length(x)
  a1 = -(n/2)*log(2*pi) - (n/2)*log(sig2)
  a2 <- -1/(2*sig2)
  y = (x-mu)^2
  
  #cat("\n n :",n)
  cat("\n mu :",mu)
  cat(" sig2 :", sig2)
  logLik = a1+a2*sum(y)
  return(logLik) 
  
}

mle_normal <- maxLik(loglik_normal, start = c(mu=1,sigma=1),x=x)
plot
summary(mle_normal)

mean(x)
var(x)*(n-1)/n
sd(x)


# varying the parameter start values

mle_normal2 <- maxLik(loglik_normal, start = c(0,1),x=x)
summary(mle_normal2)

mle_normal3 <- maxLik(loglik_normal, start = c(-1,1),x=x)    # took longer 
summary(mle_normal3)




# Grdients
logLikGrad <- function(theta,x){
  mu <- theta[1]
  sig2 <- theta[2]
  n <- length(x)
  logLikGradValues <- numeric(2)
  
  logLikGradValues[1] <- sum((x - mu)/sig2^2)
  logLikGradValues[2] <- -n/sig2 + sum((x - mu)^2/sig2^3)
  return(logLikGradValues)
} 

mle_Grad_normal <- maxLik(logLik = loglik_normal,grad =logLikGrad, start = c(1,1),x=x)
summary(mle_normal)



#########################################################################
#########################################################################

# Estimating Single parameter
# Poisson distribution

set.seed(100)
x_poisson <- rpois(n = 10000, lambda = 2)



# log likelihood function
poissonLogLik <- function(lambda,data){
  #total number of observations
  n = length(data)
  #equation
  print(lambda)
  logLik = (sum(data)*log(lambda)-n*lambda)
  #print(logLik)
  #plot(lambda,logLik)
}


mle_poisson <- maxLik(poissonLogLik,start = c(3),data=x_poisson)   # after 8 iterations
summary(mle_poisson)
mean(x_poisson)

# changing the staring values for the parameters
mle_poisson2 <- maxLik(poissonLogLik,start = c(1),data=x_poisson)  # after 9 iterations
summary(mle_poisson2)
mean(x_poisson)

data1 <- c(1,2,4,2,32,2,34,2,287)
b <- maxLik(poissonLogLik,start = c(3),data=data1)
summary(b)


lambda = 2
plot(lambda, sapply(X=lambda, FUN = function(lambda) poissonLogLik(lambda,data = x_poisson)))
  

#########################################################################
#### Poisson distribution is specified by one parameter : lambda
# This nummber equals the mean and variance
# When lambda increases to sufficiently large bumbers the normal distrbution(lambda,lambda)
# can be used to approximate the poisson distribution
#########################################################################
# Use poisson disatribution to describe hoe many times an event occurs in a finite 
# observation space

# The Poisson distribution is often used in quality control, reliability/survival studies, and insurance.
# 
# Data are counts of events(nonnegative integers with no upper bounds)
#########################################################################
# generate correct poisson random numbers 
# first 2, then 50 and 100
# plot log-likelihood




### To perform the experimennt we have to keep lambda constant
# and then plot the log-likelihood against lambda
x_poisson2 <- rpois(n = 2, lambda = 3 )
x_poisson2          # 4 0
mle_poisson2 <- maxLik(poissonLogLik,start = c(1),data=x_poisson2)   # after 1 iterations
summary(mle_poisson2)
coef(mle_poisson2)
mean(x_poisson2)


# for c = 5 : 
# Estimate : 2; std.error: 1; t-value: 2;  Pr(> t) : 0.0455     log-likelihood: -1.227411     # warnings produced
# for c = 4
#          : 2           : 0.9997    :2.001        : 0.0454                   :    # warnings produced
# for c=3      


# genberating random 50 numbers and estimating the lambda for them then takjing the monte carlo simulations
# one thing we can do to randomly change the start is by setting a random generator to select a number
# from 1-50; 

# n=50
x_poisson50 <- rpois(n = 50, lambda = 10)

x_poisson50          # 
# sum(x_poisson50) = 470
mle_poisson50 <- maxLik(poissonLogLik,start = c(5),data=x_poisson50)   # after 6 iterations
mle_poisson50 <- maxLik(poissonLogLik,start = c(8),data=x_poisson50)   # after 4 iterations
mle_poisson50 <- maxLik(poissonLogLik,start = c(9),data=x_poisson50)   # after 3 iterations        # 18 warnings
mle_poisson50 <- maxLik(poissonLogLik,start = c(11),data=x_poisson50)   # after 4 iterations
mle_poisson50 <- maxLik(poissonLogLik,start = c(13),data=x_poisson50)   # after 4 iterations       # from 13 to 7 to 9
summary(mle_poisson50)
coef(mle_poisson50)
mean(x_poisson50)


# log-likelihood: 583.1336 
# for c = 5 : 
# Estimate : 9.4000; std.error: 0.4421; t-value: 21.26;  Pr(> t) : <2e-16           # warnings produced
# for c = 8
#          :             : 0.4471               :21.02           :                   :    # warnings produced
# for c=11                                      : 22.19 




# 1793215
# store the intermediate lamnbda values and plot them with respect to value of n to check
# at what point they jumoed or divergeed to get to the lambda value (some sort of gradient descent or ascent)







