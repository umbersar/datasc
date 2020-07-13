#this is the cost func
J <- function(X, y, theta) {
  source("prediction.R")
  m = nrow(X)#length(y) or length(theta) would have been same
  (2*m)^-1 * sum((H(X,theta) - y)^2)
}


#this is a wrapper for submission system
computeCost <- function(X, y, theta) {
  J(X, y, theta) 
}

