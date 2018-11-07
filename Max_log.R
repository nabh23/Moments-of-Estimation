set.seed(100)
x <- rpois(n = 10000, lambda = 2)
#hist(x,col = 'lightblue')

getwd()

# the log-likelihood is the function of lambda and the data

poissonLogLik <- function(lambda,data){
  #total number of observations
  n = length(data)
  #equation
  print(lambda)
  logLik = (sum(data)*log(lambda)-n*lambda)
}
#########################################################
hist(x)
mle_x <- optim(par=1,fn=poissonLogLik,data=x,control=list(fnscale=-1))

# observed
table(x)
# expected
expected = dpois(x=min(x):max(x), lambda=mle_x$par)*10000


# observed
hist(x,col='lightblue')
# expected
lines(min(x):max(x),expected,col='red',lwd=2)

#########################################################
data <- 1:5
#maxLik(poissonLogLik,start = c(3),data=data)
a <- maxLik(poissonLogLik,start = c(3),data=data)
summary(a)
mean(data)

data1 <- c(1,2,4,2,32,2,34,2,2)
b <- maxLik(poissonLogLik,start = c(3),data=data1)
summary(b)
mean(data1)

c <- maxLik(poissonLogLik,start = c(3), data =x)
summary(c)
mean(x)

#observed x
table(x)

#dpois(x=1:5,lambda = 9)
# Case1 diifreent lambdas

dat <- c(0,1,3,4)
lam <- c(2,1,3,4)
p <- poissonLogLik(lam[1],data = dat)

#p1 <- poissonLogLik(lam[2],data = dat[1])
#p2 <- poissonLogLik(lam[1],data = dat[2])
#p3 <- poissonLogLik(lam[3], data = dat[3])
#p4 <- poissonLogLik(lam[4], data = dat[4])

# lambda =1 for spearate data points
#p11 <- poissonLogLik(lam[1], data = dat[1])
#p21 <- poissonLogLik(lam[1], data = dat[2])
#p31 <- poissonLogLik(lam[1], data = dat[3])
#p41 <- poissonLogLik(lam[1], data = dat[4])

#p12 <- poissonLogLik(lam[2], data = dat[1])
#p22 <- poissonLogLik(lam[2], data = dat[2])
#p32 <- poissonLogLik(lam[2], data = dat[3])
#p42 <- poissonLogLik(lam[2], data = dat[4])


#p1

#lambda <-2
p

#s <- sum(p1,p2,p3,p4)
#s

#plot(lam,c(p,p1,p2,p3,p4))
#plot(c(2,2),c(p,s))


# Case 2:

p_lam2 <- poissonLogLik(lambda = 2,data = dat)
p_lam1 <- poissonLogLik(lambda = 1, data = dat)
p_lam_onehalf <- poissonLogLik(lambda = 1.5, data = dat)
p_lam_half <- poissonLogLik(lambda = 2.5, data = dat)

p_lam2
p_lam1
p_lam_onehalf
p_lam_half

#par(mfrow=c(1,3))

plot(c(2,1,2,3,4),c(p,p1,p2,p3,p4))
plot(c(2,1,1.5,2.5),c(p_lam2,p_lam1,p_lam_onehalf,p_lam_half))
plot(lam,c(p11,p12,p13,p14))




#plot(c(2,1,1.5,2.5),c(p_lam2,p_lam1,p_lam_onehalf,p_lam_half))

p_lam2 <- poissonLogLik(lambda = 2,data = x)
p_lam1 <- poissonLogLik(lambda = 1, data = x)
p_lam_onehalf <- poissonLogLik(lambda = 1.5, data = x)
p_lam_half <- poissonLogLik(lambda = 2.5, data = x)


plot(c(2,1,1.5,2.5),c(p_lam2,p_lam1,p_lam_onehalf,p_lam_half))

################################################################
######### with non-uniform spacing in the data #################

y = c(1,20,25,30,60,65,90,150)
Y_mle <- maxLik(poissonLogLik,start = c(3),data=y)
summary(y_mle)
### warnings : In log(lambda) NaNs produced for start = 10, 3


y_mle2 = optim(par=1,fn=poissonLogLik,data=y,control=list(fnscale=-1))
# $value = 1327.235  ; par = 55.15

## calculating log likelihood for different lambdas
p_y2 <- poissonLogLik(lambda = 2,data = y)
p_y1 <- poissonLogLik(lambda = 1,data = y)
p_y3 <- poissonLogLik(lambda = 3,data = y)
p_y4 <- poissonLogLik(lambda = 25,data = y) 
p_y10 <- poissonLogLik(lambda = 10,data = y)
p_y55 <- poissonLogLik(lambda = 55,data = y)
p_y60 <- poissonLogLik(lambda = 60,data = y)

plot(c(2,1,3,25,10,55,60),c(p_y2,p_y1,p_y3,p_y4,p_y10,p_y55,p_y60))

# observed 
table(y)

expected = dpois(x=min(y):max(y), lambda=y_mle2$par)



##########################################################
############## Sample of y ###############################

sample_y <- sample(y, 100, replace=TRUE)

hist(sample_y,col='lightblue')

mle_samy <- maxLik(poissonLogLik,start = c(55.08),data=sample_y)
# 16 warnings

mle2_samy = optim(par=1,fn=poissonLogLik,data=sample_y,control=list(fnscale=-1))


p_samy5 <- poissonLogLik(lambda = 5, data = sample_y)
p_samy10 <- poissonLogLik(lambda = 10, data = sample_y)
p_samy25 <- poissonLogLik(lambda = 25, data = sample_y)
p_samy35 <- poissonLogLik(lambda = 35, data = sample_y)
p_samy50 <- poissonLogLik(lambda = 50, data = sample_y)
p_samy55 <- poissonLogLik(lambda = 55, data = sample_y)
p_samy65 <- poissonLogLik(lambda = 65,data = sample_y)

plot(c(5,10,25,35,50,55,65),c(p_samy5,p_samy10,p_samy25,p_samy35,p_samy50,p_samy55,p_samy60))

#mle2_samy = optim(par=1,fn=poissonLogLik,data=sample_y,control=list(fnscale=-1),method = "Brent",lower = 0.001,upper = 100000)
# value : 16568.39; $par: 55.07


###############################################################################
#####################Random Numbers############################################


rand_dat <-runif(10000, min = 0, max=1000)
hist(rand_dat,col='lightblue')


mle_rand <- maxLik(poissonLogLik,start = c(1),data=data_rand)
# 

mle2_rand = optim(par=1,fn=poissonLogLik,data=rand_dat,control=list(fnscale=-1),method = "Brent",lower = 0.001,upper = 100000)
## Neldar-Mead not working; using Brent

## par : 500.825;  value : 2612792



###########################################################################
###########################################################################
data_poss <-rpois(n = 10000, lambda = 55)
hist(data_poss)

mle_poss <- maxLik(poissonLogLik,start = c(54),data=data_poss)
# 36 warnings for start = 54

summary(mle_poss)
mle2_poss = optim(par=1,fn=poissonLogLik,data=data_poss,control=list(fnscale=-1))
# par: 54.94; value : 1651709




p_poss5 <- poissonLogLik(lambda = 5, data = data_poss)
p_poss10 <- poissonLogLik(lambda = 10, data = data_poss)
p_poss25 <- poissonLogLik(lambda = 25, data = data_poss)
p_poss35 <- poissonLogLik(lambda = 35, data = data_poss)
p_poss50 <- poissonLogLik(lambda = 50, data = data_poss)
p_poss55 <- poissonLogLik(lambda = 55, data = data_poss)
p_poss65 <- poissonLogLik(lambda = 65,data = data_poss)


plot(c(5,10,25,35,50,55,65),c(p_poss5,p_poss10,p_poss25,p_poss35,p_poss50,p_poss55,p_poss65))


########################################################################


hist(x,col='lightblue')

mle_x <- maxLik(poissonLogLik,start = c(3),data=x)
summary(mle_x)


p_x2 <- poissonLogLik(lambda = 2,data = x)
p_x1 <- poissonLogLik(lambda = 1,data = x)
p_x5 <- poissonLogLik(lambda = 5,data = x)
p_x25 <- poissonLogLik(lambda = 25,data = x)
p_x10 <- poissonLogLik(lambda = 10,data = x)
p_x55 <- poissonLogLik(lambda = 55,data = x)
#p_x60 <- poissonLogLik(lambda = 60,data = x)

plot(c(2,1,5,25,10,55),c(p_x2,p_x1,p_x5,p_x25,p_x10,p_x55))



################################################################

par(mfrow=c(2,2))

# custom
plot(c(2,1,3,25,10,55,60),c(p_y2,p_y1,p_y3,p_y4,p_y10,p_y55,p_y60))

# sample from custom
plot(c(5,10,25,35,50,55,65),c(p_samy5,p_samy10,p_samy25,p_samy35,p_samy50,p_samy55,p_samy60))

# poisson
plot(c(5,10,25,35,50,55,65),c(p_poss5,p_poss10,p_poss25,p_poss35,p_poss50,p_poss55,p_poss65))

plot(c(2,1,5,25,10,55),c(p_x2,p_x1,p_x5,p_x25,p_x10,p_x55))

##########################################################################

set.seed(10)
x_dat <- rpois(n = 100, lambda = 2)

p_x = poissonLogLik(lambda = 2, data = x_dat)



# simulation function:

sim_function <- function(lambda,data){

for (count in data) {
 
}
  
    
}




#hist(c(p1,p2,p3,p4))

#dat_lambda <- c(3,4,5)
#y <- rpois(n=10000, lambda)  

# 



y <- rpois(n=10000,lambda)
y <- rpois(n=10000)
dat_lambda <- c(3,4,5)
y <- rpois(10000,dat_lambda)