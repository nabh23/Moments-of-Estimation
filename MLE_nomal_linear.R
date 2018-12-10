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
  return(logLik)
}


mle_poisson <- maxLik(poissonLogLik,start = c(3),data=x_poisson)   # after 8 iterations
summary(mle_poisson)
mean(x_poisson)

# changing the staring values for the parameters
mle_poisson2 <- maxLik(poissonLogLik,start = c(1),data=x_poisson)  # after 9 iterations
summary(mle_poisson2)
mean(x_poisson)


############################################################################
################## Plotting log-likelihood #################################

#lambda = 2
#plot(lambda, sapply(X=lambda, FUN = function(lambda) poissonLogLik(lambda,data = x_poisson)))
set.seed(100)

data1 <- c(1,2,4,2,32,2,34,2,287)
b <- maxLik(poissonLogLik,start = c(3),data=data1)
summary(b)


lambda <- c(2,3,4)
lambda <- seq(0.5,5,length.out = 100)
#plot(lambda, sapply(X=lambda, FUN = function(lambda) poissonLogLik(lambda,data = x_poisson)))

vllh <- Vectorize(poissonLogLik,"lambda")
vllh(c(2,3,4),x_poisson)
plot(lambda, vllh(lambda,x_poisson), type="l")


# Generating 2 numbers
x_poisson2 <- rpois(n = 2, lambda = 3 )     # Here the dataset is created with lambda=3
# But my estimate is 2; Also the log-likelihood is maximum for lambda = 2
x_poisson2
summary(maxLik(poissonLogLik,start = c(3),data=x_poisson2))
lambda2 <- c(2,3,4)
#vllh2 <- Vectorize(poissonLogLik,"lambda2")
vllh(c(2,3,4), x_poisson2)
# plot of log-lokelihood for 2 random numbers with different lambdas
plot(lambda, vllh(lambda,x_poisson2), type="l")



# Generating 50 numbers

# 1. With lambda = 3


x_poisson50_1 <- rpois(n = 50, lambda = 3)
# using the same range of lambdas as before to estimate
# i.e c(2,3,4)
mean(x_poisson50_1)
vllh(c(2,3,4), x_poisson50_1)
plot(lambda, vllh(lambda,x_poisson50_1), type="l")

# 2. With lambda = 5
x_poisson50_2 <- rpois(n = 50, lambda = 5)
x_poisson50_2
mean(x_poisson50_2)
summary(maxLik(poissonLogLik,start = c(3),data=x_poisson50_2))

lambda3 <- c(3,4,5,7)
vllh(c(3,4,5,7), x_poisson50_2)
plot(lambda3, vllh(lambda3,x_poisson50_2), type="l")

## Why the log-likelihoods are positive ?? Does the value of kambda during generation of
# numbers have any affect on the log-likelihood?

# 3. With lambda = 2
x_poisson50_3 <- rpois(n = 50, lambda = 2)
x_poisson50_3
mean(x_poisson50_3)
summary((maxLik(poissonLogLik,start = c(3),data=x_poisson50_3)))

lambda4 <- c(1,3,2,4,5)
vllh(c(1,3,2,4,5), x_poisson50_3)
plot(lambda4, vllh(lambda4,x_poisson50_3), type="l")


# Generating 100 random poisson distributions

# 1. With lambda = 2
x_poisson100_1 <- rpois(n = 100, lambda = 2)
x_poisson100_1
mean(x_poisson100_1)
summary((maxLik(poissonLogLik,start = c(4),data=x_poisson100_1)))

lambda5 <- c(1,2,3,4,5)
vllh(c(1,2,3,4,5), x_poisson100_1)
plot(lambda, vllh(lambda,x_poisson100_1), type="l")


# 2. with lambda = 4
x_poisson100_2 <- rpois(n = 100, lambda = 5)
x_poisson100_2
mean(x_poisson100_2)
mle_poisson100_2 <- (maxLik(poissonLogLik,start = c(4),data=x_poisson100_2))
summary(mle_poisson100_2)
#summary((maxLik(poissonLogLik,start = c(11),data=rpois(n = 100, lambda = 9))))

lambda6 <- c(2,3,5,6,7)
vllh(c(2,3,5,6,7), x_poisson100_2)
plot(lambda6, vllh(lambda6,x_poisson100_2), type="l")
#abline(v=mean(x_poisson100_2), lty=2)
abline(v=mle_poisson100_2$estimate, lty=2)
mtext(expression(hat(y)), side=1, at=(mle_poisson100_2$estimate), padj=1, cex=1.1)

# 3. with lambda = 9
x_poisson100_3 <- rpois(n = 100, lambda = 9)
x_poisson100_3
mean(x_poisson100_3)
summary((maxLik(poissonLogLik,start = c(6),data=x_poisson100_3)))
#summary((maxLik(poissonLogLik,start = c(11),data=rpois(n = 100, lambda = 9))))

lambda7 <- c(5,7,8,9,11,12)
vllh(c(5,7,8,9,11,12), x_poisson100_3)
plot(lambda7, vllh(lambda7,x_poisson100_3), type="l")
#abline(v=mean(x_poisson100_3), lty=2)
abline(v=mle_poisson100_3$estimate, lty=2)
mtext(expression(hat(y)), side=1, at=mean(x_poisson100_3$estimate), padj=1, cex=1.1)


########################################################################
## Now repeating same steps for random non-uniformly spaced numbers

# Function to generate random non-uniformly spaced numbers

mysample <- function(n, lower, upper, space) {
  b <- ceiling((upper - lower + 1) / (space - 1))
  bs <- sample(seq(2, b - 1, by = 2), n - 1)
  gr <- split(setdiff(1:b,bs), cumsum(c(0, diff(setdiff(1:b, bs))) != 1))
  out <- sapply(gr, function(x) (x[1] - 1) * (space - 1) + ceiling(runif(1) * length(x) * (space - 1)))
  out[n] <- min(out[n], upper)
  out
}



# 1. generating 2 numbers
#set.seed(123)
x_random2 <- replicate(min(diff(mysample(100, 1, 100000, 5))), n = 2)
x_random2
min(x_random2)
mean(x_random2)
mle_ran_2 <- (maxLik(poissonLogLik,start = c(52),data=x_random2))
summary(mle_ran_2)
#summary((maxLik(poissonLogLik,start = c(11),data=rpois(n = 100, lambda = 9))))

lambda_ran_1 <- seq(50,54,.5)
vllh(seq(50,54, .5), x_random2)
plot(lambda_ran_1, vllh(lambda_ran_1,x_random2), type="l", xlab="lambda", ylab="log-likelihood")
#abline(v=mean(x_random2), lty=2)
abline(v=mle_ran_2$estimate, lty=2)
mtext(expression(hat(y)), side=1, at=mean(mle_ran_2$estimate), padj=1, cex=1.1)


# 2. generating 50 numbers

x_random50 <- replicate(min(diff(mysample(100, 1, 100000, 5))), n = 50)
x_random50
min(x_random50)
mean(x_random50)
mle_ran_2 <- (maxLik(poissonLogLik,start = c(34),data=x_random50))
summary(mle_ran_2)
#summary((maxLik(poissonLogLik,start = c(11),data=rpois(n = 100, lambda = 9))))

lambda_ran_2 <- seq(29.5,36,.5)
vllh(seq(29.5,33,.5), x_random50)
plot(lambda_ran_2, vllh(lambda_ran_2,x_random50), type="l",xlab="lambda", ylab="log-likelihood")
#abline(v=mean(x_random2), lty=2)
abline(v=mle_ran_2$estimate, lty=2)
abline(h=mle_ran_2$maximum, lty=1)
mtext(expression(hat(y)), side=1, at=mean(mle_ran_2$estimate), padj=1, cex=1.1)


# 3. generating 100 numbers

x_random100 <- replicate(min(diff(mysample(100, 1, 100000, 2))), n = 100)
x_random100
min(x_random100)
mean(x_random100)
mle_ran_3 <- (maxLik(poissonLogLik,start = c(36),data=x_random100))
summary(mle_ran_3)
#summary((maxLik(poissonLogLik,start = c(11),data=rpois(n = 100, lambda = 9))))

lambda_ran_3 <- seq(29,36,.6)
vllh(seq(31,36,.6), x_random100)
plot(lambda_ran_3, vllh(lambda_ran_3,x_random100), type="l",xlab="lambda", ylab="log-likelihood")
#abline(v=mean(x_random2), lty=2)
abline(v=mle_ran_3$estimate, lty=2)
abline(h=mle_ran_3$maximum, lty=1)
mtext(expression(hat(y)), side=1, at=mean(mle_ran_3$estimate), padj=1, cex=1.1)



#####################
x_random200 <- replicate((diff(mysample(100, 1, 10000, 2.5))), n = 1)
x_random200
min(x_random200)
mean(x_random200)
mle_ran_4 <- (maxLik(poissonLogLik,start = c(98),data=x_random200))
summary(mle_ran_4)



lambda_ran_4 <- seq(97,103,.6)
vllh(seq(97,102,.6), x_random200)
plot(lambda_ran_4, vllh(lambda_ran_4,x_random200), type="l",xlab="lambda", ylab="log-likelihood")
#abline(v=mean(x_random2), lty=2)
abline(v=mle_ran_4$estimate, lty=2)
abline(h=mle_ran_4$maximum, lty=1)
mtext(expression(hat(y)), side=1, at=mean(mle_ran_4$estimate), padj=1, cex=1.1)

########### The maxLik or MLE is able to estimate data with discrepency in numbers
# 


## We sampled 3000 with replacement from a set of 586







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
hist(x_poisson50)
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


# 100 random numbers
set.seed(100)
x_poisson100 <- rpois(n = 100, lambda = 10)

x_poisson100          # 
hist(x_poisson100)
mean(x_poisson100)
# sum(x_poisson50) = 470
mle_poisson100 <- maxLik(poissonLogLik,start = c(5),data=x_poisson100)   # after 5 iterations
mle_poisson100 <- maxLik(poissonLogLik,start = c(8),data=x_poisson100)   # after 4 iterations
mle_poisson100 <- maxLik(poissonLogLik,start = c(9),data=x_poisson100)   # after 3 iterations        # 18 warnings
mle_poisson100 <- maxLik(poissonLogLik,start = c(11),data=x_poisson100)   # after 4 iterations
mle_poisson100 <- maxLik(poissonLogLik,start = c(13),data=x_poisson100)   # after 5 iterations       # from 13 to 7 to 9
summary(mle_poisson100)
coef(mle_poisson100)
mean(x_poisson100)


# Estimate : 1256.735;      
# estimate : 9,8000           std.erros: 0.3198     tvalue: 30.64       Pr(> t): <2e-16
#                                         0.3162                   31
#                                        0.2996             32.71
                                                            

x_poisson_random <- c(12,5,12,9,8)
hist(x_poisson_random)
poissonLogLik(10, x_poisson_random)
sum(x_poisson_random)
mean(x_poisson_random)
maxLik(poissonLogLik, start = c(9), data = x_poisson_random)



# Now, generating non-uniform numbers , i.e, numbers with gaps in between

