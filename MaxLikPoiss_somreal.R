set.seed(100)
x <- rpois(n = 10000, lambda = 2)
hist(x,col = 'lightblue')

getwd()

# the log-likelihood is the function of lambda and the data

poissonLogLik <- function(lambda,data){
  #total number of observations
  n = length(data)
  #equation
  logLik = (sum(data)*log(lambda)-n*lambda)
}




# using optim fn to find the max

mle <- optim(par=1,fn=poissonLogLik,data=x,control = list(fnscale=-1))
mle

##Error## for optim
# Warning message:
#In optim(par = 1, fn = poissonLogLik, data = x, control = list(fnscale = -1)) :
#one-dimensional optimization by Nelder-Mead is unreliable:
#use "Brent" or optimize() directly

# observed value
table(x)

# expected values
expected = dpois(x= min(x):max(x),lambda = mle$par)*10000   # x = 0:9, par$mle = 2.00019
expected

#observed
hist(x,col = 'lightblue')

#expected
lines(min(x):max(x), expected, col= 'red', lwd=2)
