plague <- c(20, 40, 31, 17, 10, 12, 14, 4, 10, 12, 9, 16, 8, 5, 4, 4, 9)
plot(1982:1998,plague,
     type="l",
     xlab="Year",
     ylab="Plague Cases")
f <- function(x,lambda) -sum(log(dpois(x,lambda)))
aaa <- optimize(f, c(1,20), x=plague)
aaa$min


ll <- rep(0,20)

for(i in 1:20) ll[i] <- -sum(log(dpois(x_poisson,i)))
plot(1:20, ll,type="l", log="y" ,xlab="year", ylab="negative log-likelihood")
abline(h=aaa$objective,lty=2)
abline(v=aaa$min, lty=2)
mtext(expression(hat(y)), side=1, at=aaa$min, padj=1, cex=1.1)

##########################################

set.seed(100)
x2 <- rpois(n = 2, lambda = 3 )
x50<- rpois(n=50, lambda =3)

input <- c(x2,x50)

ll <- rep(0,2)

for(i in input) ll[i] <- -sum(log(dpois(i,lambda)))
plot(1:8, ll, type = "l")
