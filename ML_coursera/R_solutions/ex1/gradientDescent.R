gradientDescent <- function(X, y, theta, alpha, iterations) {
  source("prediction.R")#includes the hypothesis function
  m = nrow(X)#length(y) or length(theta) would have been same
  
  #save theta history as well as cost history for that theta so that we can plot the route on a countour map
  j_history <- vector(mode="numeric", length=iterations+1)#j_history = rep(0,iterations)
  theta_history <- matrix(data=0,nrow=iterations+1,ncol=length(theta))#initilize a matrix with 0's
  
  #store the initial cost and theta in history var
  j_history[1] = J(X,y, theta)
  theta_history[1,] = theta
  
  for (i in (2:iterations+1)) {
    slope = m^-1 * (t(X) %*% (H(X,theta)-y))#the derivation term
    theta = theta - (alpha)*slope
    
    j_history[i] = J(X,y, theta)
    theta_history[i,] = theta
  }
  
  #theta_history[iterations+1] is the last theta at the end of iterations
  list(theta = theta_history[iterations+1,], j_history=j_history, theta_history=theta_history)#this will be returned
  # ------------------------------------------------------------
}
