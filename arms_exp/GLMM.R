rm(list=ls())
dataDir = "/media/sf_Ubuntu/Rcpptest/"
setwd(dataDir)

library(Rcpp)
library(RcppArmadillo)
library(HI)


# read in data
#---------------------------------------------------
dat <- read.table("./mydata.txt", header=T)
I <- length(unique(dat$ID))
myID <- unique(dat$ID)
nobs <- NULL # number of observations from each subject
subject <- NULL # subject number of each observation
for (i in 1:I) {
  temp <- subset(dat, ID == myID[i])
  temp.nrow <- nrow(temp)
  nobs <- c(nobs, temp.nrow)
  subject <- c(subject, rep(i,temp.nrow))
}
ncycles <- 5000
p <- 2
design.matrix <- cbind(rep(1,I), dat$x)
y <- dat$y
myshape <- 0.5*I-1 # the alpha parameter for the full conditional of parameter tau
#---------------------------------------------------

temp.beta <- rep(1, p); temp.u <- rep(0.1, I); temp.tau <- 1

sourceCpp('arms.cpp')
#--------------------------------------------------------------------------
time <- Sys.time()
out = sampling(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles)
beta1 = out$beta
tau1 = out$tau
u1 = out$u
time1 = Sys.time()-time
#------------------------------------------------

time <- Sys.time()
logden.u <- function(x, mybeta, myy, mycovariate, mytau) {
  temp <- mycovariate %*% as.matrix(mybeta, ncol=1) + x
  return(sum(myy * temp - log(1+exp(temp))) - 0.5*x^2*mytau)
}

logden.beta <- function(x, myu) {
  myu.long <- as.vector(mapply(rep, myu, nobs))
  temp <- as.vector(design.matrix %*% as.matrix(x, ncol=1)) + myu.long
  return(sum(y*temp - log(1+exp(temp))))
}


beta2 <- NULL
tau2 <- NULL
u2 <- NULL
for (n in 1:ncycles) {
  temp.beta <- arms(temp.beta, logden.beta, function(x,...) all(x>(-10))*all(x<(10)), 1, myu=temp.u)
  temp.tau <- rgamma(1, shape=myshape, rate=0.5*sum(temp.u^2))
  beta2 = cbind(beta2, temp.beta)
  tau2 = c(tau2, temp.tau)
  
  # update random intercept U_i
  save.u <- NULL; counter <- 0
  for (i in 1:I) {
    start <- counter + 1
    end <- counter + nobs[i]
    temp <- arms(temp.u[i], logden.u, function(x,...) (x>(-100))*(x<100), 1,
                 mybeta=temp.beta, myy=y[start:end], mycovariate=design.matrix[start:end,],
                 mytau=temp.tau)
    save.u <- c(save.u, temp)
    counter <- end
  }
  temp.u <- save.u; rm(save.u)
  u2 = cbind(u2,temp.u)
}
time2 = Sys.time()-time

#---------------------------------------------------
print(time1)
print(time2)

# > print(time1)
# Time difference of 1.177525 mins
# > print(time2)
# Time difference of 36.14044 mins

par(mfrow=c(1,2))
plot(density(beta1[1,-c(1:400)]), xlab = "", main = "beta0 cpp")
plot(density(beta2[1,-c(1:400)]), xlab = "", main = "beta0 R")
plot(density(beta1[2,-c(1:400)]), xlab = "", main = "beta1 cpp")
plot(density(beta2[2,-c(1:400)]), xlab = "", main = "beta1 R")

plot(density(tau1[-c(1:400)]), xlab = "", main = "tau cpp")
plot(density(tau2[-c(1:400)]), xlab = "", main = "tau R")

u1.mean = apply(u1[,-c(1:400)],1,mean)
u2.mean = apply(u2[,-c(1:400)],1,mean)
par(mfrow=c(1,2))
plot(density(u1.mean), xlab = "", main = "us cpp")
plot(density(u2.mean), xlab = "", main = "us R")

save(beta1,beta2, tau1,tau2, u1,u2, file="samples.rdata")