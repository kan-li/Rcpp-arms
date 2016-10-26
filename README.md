# Rcpp-arms

The <a href="https://cran.r-project.org/web/packages/Rcpp/index.html" target="_blank">Rcpp package</a> provides R functions and a C++ library facilitating the integration of R and C++.

The code here provides an example of MCMC for generalized linear mixed model using Adaptive Metropolis rejection sampling (AMRS). 

The test folder include a test example for Rcpp, RcppArmadillo and parallelization in Rcpp via OpenMP. Three ways are used to calculate the inner product of a vector: R operation '%*%', user defined function 'inner' coded in Rcpp, and parallelize the 'inner' function with multiple cores. A small set of micro-benchmarks in a variety of methods are conducted.

The arms-exp folder include the example of MCMC for generalized linear mixed model using Adaptive Metropolis rejection sampling (AMRS). Two methods are compared. First  is using the arms function in R Hi package. Second is writing log density functions and sampling process all in Rcpp/C. The pure cpp code is 30 times faster than R code. 

The arms_openMP_exp folder include the same example but parallelized the sampling process for random effect. Usually this part is the most time consuming one when the number of subject is large. However, sampling random effect for each subject is independent can be parallelized. However, I didn't find a big time saving when the using multiple cores in this case. This may because the time spent on transferring/aggregating data is longer than the time saved for sampling. Further testing may needed.

If you are new to Rcpp then the following resources provide a helpful introduction:

- Rcpp: Seamless R and C++ Integration
- [Rcpp website](http://www.rcpp.org/)
- [Rcpp book](http://www.rcpp.org/book/)
- [Rcpp Quick Reference Guide](http://cran.rstudio.com/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf)
- [High performance functions with Rcpp](http://adv-r.had.co.nz/Rcpp.html)
- [Introduction to Rcpp Attributes](http://cran.rstudio.com/web/packages/Rcpp/vignettes/Rcpp-attributes.pdf)
- [Gallery of examples](http://gallery.rcpp.org/)
- [Armadillo C++ linear algebra library](http://arma.sourceforge.net/docs.html)