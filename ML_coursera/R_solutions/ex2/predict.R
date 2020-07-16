
predict <- function(theta, X){
  #returns 0 or 1 using learned logistic regression
  prediction_before_logistic_transformation <- X%*%theta
  
  #call the sigmoid hypothesis func
  logistic_regression_hypothesis_val = h(prediction_before_logistic_transformation)
  #.5 is the threshold
  ifelse(logistic_regression_hypothesis_val>=.5, 1, 0)
}

