library(HI)

# A modified version of arms with as.integer(n.sample) = 1 

my.arms = function (y.start, myldens, indFunc, ...) {   
  #bounds <- y.start + convex.bounds(y.start, dir = 1, indFunc, ...)
  bounds <- c(-100,100)
  if (diff(bounds) < 1e-07) 
    y.sample <- rep(y.start, 1)
  else {
    f <- function(x) myldens(x, ...)
    y.sample <- .Call("arms", bounds, f, y.start, as.integer(1), new.env())
  }
  return(y.sample)
}