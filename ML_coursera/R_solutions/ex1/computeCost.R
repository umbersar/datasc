#this is the cost func
J <- function(X, y, theta) {
  m = nrow(X)#length(y) or length(theta) would have been same
  (2*m)^-1 * sum((H(X,theta) - y)^2)
}
