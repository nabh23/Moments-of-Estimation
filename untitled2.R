nobs=1000

x1 = rnorm(nobs)*.15^.5
x2 = rnorm(nobs)*.35^.5
u = rnorm(nobs)*(.5)^.5

inside = 2*x1+2*x2+u

# Using basic
# Var(inside) = (.15*2)^2*var(N(0,1)) + (.35*2)^2*var(N(0,1)) + (.5^.5)^2*var(N(0,1)) = .15+.35+.5 =1


# The probability of have a success is equal the cumulative density function at the inside value.
p = pnorm(inside)

y = rbinom(nobs,1,p)

mydata = data.frame(y=y,x1=x1,x2=x2)


# The estimates we would like to duplicate are:
summary(glm(y~x1+x2, family=binomial(link="probit"), data=mydata))


# function that takes the value and maximizes over

myprobit <- function(theta) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  
  log.lik <- sum(log(pnorm((-1)^(1-mydata$y)*(beta0 + beta1*mydata$x1 + beta2*mydata$x2))))
  
  return(log.lik)
}

myprobit(c(1,1,1))


# Let's see what the probit log likelihoods look like as we vary B1 from -1 to 2
b.range <- seq(-1,2.5,.1)
b0.loglik <- b1.loglik <- b2.loglik <- 0*b.range
for (i in 1:length(b.range)) {
  b0.loglik[i] <- myprobit(c(b.range[i],1,1))
  b1.loglik[i] <- myprobit(c(1,b.range[i],1))
  b2.loglik[i] <- myprobit(c(1,1,b.range[i]))
}
  

# I will define a quick min max function to help setting up the range of plots
minmax <- function(x) c(min(x),max(x))

plot(minmax(b.range), minmax(c(b0.loglik, b1.loglik, b2.loglik)),
     type="n", main="Log-Likelihood Around the point (1,1,1)",
     ylab = "log likelihoods", xlab= ~ beta )

lines(c(1,1),minmax(c(b0.loglik, b1.loglik, b2.loglik)), col="green", lwd=2)

lines(b.range, b0.loglik, col="black", lwd=2)
lines(b.range, b1.loglik, col="blue", lwd=2)
lines(b.range, b2.loglik, col="purple", lwd=2)


optim(c(1,1,1), myprobit, control=list(fnscale = -1))

