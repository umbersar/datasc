J  <- function(X, y, theta) {
  m = nrow(X)#number of rows in training dataset
  #(2*m)^-1 * sum((H - y)^2)#for reference, this was MSE cost func for Linear regression
  
  cost <- -1*( (1/m) * ( t(y) %*% log(H(X %*% theta)) + t(1-y) %*% log(1-H(X %*% theta)) ) )
}

#The gradient/slope of the cost func w.r.t. to the parameters thetas.
grad <- function(X, y, theta) {
  grad <-  m^-1 * (t(X) %*% (H(X %*% theta)-y))
  grad
}

#to be callable by optim function, you have to modify the cost func J to be written to be written
#as closure on X and Y and callable using theta:
J_for_optim  <- function(X, y) {
  function(theta) {
    m <- length(y) # number of training examples
    J <- -1*( (1/m) * ( t(y) %*% log(H(X %*% theta)) + t(1-y) %*% log(1-H(X %*% theta)) ) )
  }
}

#The gradient/slope of the cost func w.r.t. to the parameters thetas. This func has to be
#to be written to be rewritten as closure on X and Y and callable using theta so that optim
#can use it:
grad_for_optim <- function(X, y) {
  function(theta) {
    grad <-  m^-1 * (t(X) %*% (H(X %*% theta)-y))
    grad
  }
}

