# Rcpp-arms

The <a href="https://cran.r-project.org/web/packages/Rcpp/index.html" target="_blank">Rcpp package</a> provides R functions and a C++ library facilitating the integration of R and C++.

The code here provides an example of MCMC for generalized linear mixed model using Adaptive rejection Metropolis sampling (ARMS). 

The test folder include a test example for Rcpp, RcppArmadillo and parallelization in Rcpp via OpenMP. Three ways are used to calculate the inner product of a vector: R operation '%*%', user defined function 'inner' coded in Rcpp, and parallelize the 'inner' function with multiple cores. A small set of micro-benchmarks in a variety of methods are conducted.

The arms-exp and arms-hybrid-exp folder include the example of MCMC for generalized linear mixed model using ARMS.Three methods are compared. The first method is pure R code, using the arms function in R Hi package. The second method is writing log density functions in C++, but sampling process is still coded in R with arms function.  The third method is coding both log density functions and sampling process in Rcpp/C++. It calls the original C code of ARMS written by  Wally Gilks. One example of how to call this function is provided on the <a href="https://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html" target="_blank">webpage</a>.  The pure R code takes 36.14 mins. The hybrid method takes 9.80 mins, and the C++ code used only 1.17 mins.

The arms_openMP_exp folder include the same example but parallelized the sampling process for random effect. Usually this part is the most time consuming one when the number of subject is large. However, sampling random effect for each subject is independent and can be parallelized. However, I didn't find a big time saving when increased the number of cores more than 3 in this case. This may because the time spent on transferring/aggregating data is longer than the time saved for sampling. Further testing may needed.

If you are new to Rcpp then the following resources provide a helpful introduction:

- Rcpp: Seamless R and C++ Integration
- [Rcpp website](http://www.rcpp.org/)
- [Rcpp book](http://www.rcpp.org/book/)
- [Rcpp Quick Reference Guide](http://cran.rstudio.com/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf)
- [High performance functions with Rcpp](http://adv-r.had.co.nz/Rcpp.html)
- [Introduction to Rcpp Attributes](http://cran.rstudio.com/web/packages/Rcpp/vignettes/Rcpp-attributes.pdf)
- [Gallery of examples](http://gallery.rcpp.org/)
- [Armadillo C++ linear algebra library](http://arma.sourceforge.net/docs.html)