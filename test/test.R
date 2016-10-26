library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

sourceCpp('test.cpp')

bench = benchmark(vecA%*%vecA,
          inner(vecA,vecA),
          sumsq_parallel(vecA,1),
          sumsq_parallel(vecA,2),
          sumsq_parallel(vecA,3),
          sumsq_parallel(vecA,4), order = NULL, replications=1000)

print(bench[,1:4])

# test replications elapsed relative
# 1           vecA %*% vecA         1000   0.721    4.682
# 2       inner(vecA, vecA)         1000   0.229    1.487
# 3 sumsq_parallel(vecA, 1)         1000   0.319    2.071
# 4 sumsq_parallel(vecA, 2)         1000   0.203    1.318
# 5 sumsq_parallel(vecA, 3)         1000   0.154    1.000
# 6 sumsq_parallel(vecA, 4)         1000   0.157    1.019