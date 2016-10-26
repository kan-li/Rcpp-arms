rm(list=ls())
dataDir = "/home/sph/kli/Rcpp/"
#setwd(dataDir)

library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")



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


sourceCpp(paste(c(dataDir,'arms_par.cpp'), collapse = ""))

brench = benchmark(sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 1),
                   sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 3),
                   sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 5),
                   sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 7),
                   sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 9),
                   sampling_par(nobs, design.matrix, y, myshape,temp.beta, temp.u, temp.tau, ncycles, 11),
                   order=NULL, replications = 50)
save(brench, file="benchmark.rdata")
print(brench[,1:4])

