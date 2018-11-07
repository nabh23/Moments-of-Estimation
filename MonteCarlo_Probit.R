# how assymptotic variance works in probit models

Probit_LL <- function(y,x,par) {
  Phi = pnorm(x %*% par)
  phi = dnorm(x %*% par)
  
  n = length(y) 
  k = length(par)
  
  # Computing the log-likelihood
  f = sum(y*log(Phi)) + sum((1-y)*log(1-Phi))
  f = -f
  
  return(f)
}


# function for th egradient

Probit_LL_g <- function (y,x,par) {
  Phi = pnorm(x %*% par) # Phi is Cumulative probability
  phi = dnorm(x %*% par) # phi is Probability Density
  
  n = length(y)           # sample size
  k = length(par)         # number of coefficients
  
  g = t(matrix(rep(phi/Phi,k),nrow=n)*x) %*% y - 
    t(matrix(rep(phi/(1-Phi),k),nrow=n)*x) %*% (1-y)
  g = -g
  
  return(g)
} 


# Monte Carlo Simulation
# assymptotic concept


library(numDeriv)           # loads the functions grad and hessian which numeriacally evaluate the gradient and hessian
#rm(list=ls())       # Clear workspace

set.seed(1)         # Set seed for random number generator

n = 100         # Set the sample size

theta0 = c(1,1,1)   # True parameter value

k = length(theta0) 

# Data generating process
x = cbind(matrix(1,n,1),matrix(rnorm(2*n,0,1),ncol=2))  # regressors

u = rnorm(n,0,1)                    # error term

y_star = x %*% theta0 + u               # latent "utility"

y = ceiling(y_star/max(abs(y_star)))    # observed outcome (y=1 if y_star >0, otherwise y=0)

dat = data.frame(x,y)




# Checkin the gradient was loaded properly or not

Probit_LL_g(y,x,theta0)

grad(function(u) Probit_LL(y,x,u),theta0)

# using optim to estimate the paraemeters

result <- optim(par = theta0, Probit_LL, y = y, x = x, gr = Probit_LL_g, 
                method = c("BFGS"), control = list(reltol=1e-9), hessian=TRUE)



