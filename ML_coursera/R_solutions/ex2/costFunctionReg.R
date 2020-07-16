J  <- function(X, y, theta, lambda) {
  m = nrow(X)#number of rows in training dataset
  #(2*m)^-1 * sum((H - y)^2)#just for reference this was MSE cost func for Linear regression
  
  #regularized cost func with added penalty on thetas (except theta0)
  cost <- -1*( (1/m) * ( t(y) %*% log(H(X %*% theta)) + t(1-y) %*% log(1-H(X %*% theta)) ) )
          + (lambda/(2*m)) * sum(theta[-1]^2)
  
  #a vector to hold values of slopes/gradients for all thetas/params
  slope = vector(mode = "numeric", length = length(theta))#could also use rep(0,length(theta))
  
  #For regularized logistic regression, we have to separate the vectorized slope/gradient 
  #equation into 2 as we only add the penalty to thetas/params other than theta0
  slope[1] = (1/m)* t(X[,1]) %*% (H(X %*% theta) - y)
  slope[-1] = (1/m)*  t(X[,-1]) %*% (H(X %*% theta) - y) + (lambda/m)*theta[-1]
  
  list(J=cost,Slope=slope)
  
}


#just a wrapper to satisfy grader
costFunction <- function(theta, X, y, lambda){
  cost_obj <- J(X,y,theta)
  list(J=cost_obj$J, grad=cost_obj$Slope)
}
