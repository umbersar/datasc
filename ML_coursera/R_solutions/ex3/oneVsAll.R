oneVsAll <- function(X, y, num_labels, lambda) {
  #   all_theta <- ONEVSALL(X, y, num_labels, lambda) trains num_labels
  #   logisitc regression classifiers and returns each of these classifiers
  #   in a matrix all_theta, where the i-th row of all_theta corresponds
  #   to the classifier for label i
  
  # Some useful variables
  m <- nrow(X)
  n <- ncol(X)
  
  #this will contain the theta/params for each of the oneVSall classifiers. Other way to put it is
  #that each row of the theta matrix will contain the theta for predicting a label
  all_theta <- matrix(0, nrow = num_labels, ncol = n + 1)
  
  # add the dummy feature to X data matrix
  X <- cbind(rep(1,m),X)
  

  #To train oneVsall classifier, we run the optimization for each of the labels. And for each label,
  #we set the training label to 0 for other labels(so it becomes a binary classification wrt that label)
  #in the training label vector
  #note here that the the number 1 is represented by 1 in training label,....till 9 and then 10 
  #represents(or is the label of) digit 0
  for (i in 1:num_labels) {
    # Set Initial theta
    initial_theta <- rep(0,n + 1)
    
    # Run optim to obtain the optimal theta. This function will return theta and the cost
    opt <- optim(par = initial_theta, fn = J_regulalized_for_optim(X,y == i,lambda),
      gr = grad_regulalized_for_optim(X,y == i,lambda),method = "BFGS"
      ,control = list(maxit = 50))
    
    theta <- opt$par
    J <- opt$value# we don't have any use for cost here but still
    cat(sprintf("Iteration %d | Min Cost: %f\n",i,opt$value))
    
    all_theta[i,] <- theta
  }
  
  all_theta
}
