#the argument list would have been function(X, y, theta, lambda)
#but since we are going to use optim to find the theta/params (by minimizing the cost), we need
#to make a closure of this func on arguments except theta/params
J_regulalized_for_optim <- function(X, y, lambda) {
  #lrCostFunction Compute cost for logistic regression with
  #regularization
  
  function(theta) {
    m = nrow(X)#number of rows in training dataset
    #(2*m)^-1 * sum((H - y)^2)#just for reference this was MSE cost func for Linear regression
    
    #regularized cost func with added penalty on thetas (except theta0)
    cost <- -1*( (1/m) * ( t(y) %*% log(H(X %*% theta)) + t(1-y) %*% log(1-H(X %*% theta)) ) )
    + (lambda/(2*m)) * sum(theta[-1]^2)
    cost
  }
}

grad_regulalized_for_optim <- function(X, y, lambda) {
  function(theta) {
    #For regularized logistic regression, we have to separate the vectorized slope/gradient 
    #equation into 2 as we only add the penalty to thetas/params other than theta0
    m = nrow(X)#number of rows in training dataset
    
    #a vector to hold values of slopes/gradients for all thetas/params
    grad = vector(mode = "numeric", length = length(theta))#could also use rep(0,length(theta))
    
    grad[1] = (1/m)* t(X[,1]) %*% (H(X %*% theta) - y)
    grad[-1] = (1/m)*  t(X[,-1]) %*% (H(X %*% theta) - y) + (lambda/m)*theta[-1]
    grad    
    
  }
}
