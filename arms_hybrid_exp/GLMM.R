rm(list=ls())

library(Rcpp)
library(RcppArmadillo)

source("myarms.R")

# library(foreach)
# library(doParallel)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)


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

sourceCpp('GLMMlden.cpp')

time <- Sys.time()

beta <- NULL
tau <- NULL
u <- NULL

start <- rep(0,I)
end <- rep(0,I)
counter <- 0
for(i in 1:I){
  start[i] <- counter + 1
  end[i] <- counter + nobs[i]
  counter=end[i];
}


for (n in 1:ncycles) {
  myu.long <- as.vector(mapply(rep, temp.u, nobs))
  
  temp.beta <- if(length(temp.beta>1)){
    y.start = temp.beta
    myldens = logden_beta
    indFunc = function(x,...) all(x>(-10))*all(x<(10))
    dim <- length(temp.beta)
    y.sample <- rbind(y.start, matrix(0, 1, dim))
    dir <- rnorm(dim)
    #bounds <- convex.bounds(y.sample[1, ], dir, indFunc, myu_long = myu.long, design_matrix=design.matrix, y=y)
    bounds = c(-10,10)
    if (diff(bounds) < 1e-07) 
      y.sample[2, ] <- y.sample[1, ]
    else {
      f <- function(x) myldens(y.sample[1, ] + x *dir, myu_long = myu.long, design_matrix=design.matrix, y=y)
      y.sample[2, ] <- y.sample[1, ] + dir * .Call("arms",bounds, f, 0, as.integer(1), new.env())
    }
    y.sample <- y.sample[-1, ]
  }
  
  temp.tau <- rgamma(1, shape=myshape, rate=0.5*sum(temp.u^2))
  beta = cbind(beta, temp.beta)
  tau = c(tau, temp.tau)
  
  # update random intercept U_i
  temp.u <- rep(0, I); counter <- 0
  for (i in 1:I) {
    temp.u[i] <- my.arms(temp.u[i], logden_u, function(x,...) (x>(-100))*(x<100), 
                 mybeta=temp.beta, myy=y[start[i]:end[i]], mycovariate=design.matrix[start[i]:end[i],],mytau=temp.tau)
  }
  

  
#   temp.u <- foreach(i=1:I, .combine = cbind) %dopar% {
#          a = my.arms(temp.u[i], logden_u, function(x,...) (x>(-100))*(x<100), 
#                      mybeta=temp.beta, myy=y[start[i]:end[i]], mycovariate=design.matrix[start[i]:end[i],],mytau=temp.tau)
#          a
#   }
  
  u = cbind(u,temp.u)
}
time = Sys.time()-time

print(time)


# Time difference of 9.796264 mins
