predictOneVsAll  <- function(all_theta, X) {
  #   p <- PREDICT(all_theta, X) computes the predictions for X using a
  #   threshold at 0.5 (i.e., if sigmoid(t(all_theta) %*% x) >= 0.5, predict 1)
  
  m <- nrow(X)#X contains the unlabeled data for which we have to predict.
  num_labels <- nrow(all_theta)#all_theta has the thetas/params for all the oneVSall classifiers
  
  # You need to return the following variables correctly
  p <- rep(0,nrow(X))#p is the prediction vector of length of m and has the prediction 
  #for every row in the unlabeled data. 
  
  # add the dummy feature to the matrix of features(X...unlabeled data)
  X <- cbind(rep(1,m), X)
  
  #apply every theta/params (each row from all_theta) to every each and every row in X
  #For, e.g., apply all the different theta sets to the first row from X. It will give you the 
  #probability of that row from X being one of the label(the thetas from all_thetas are for correct
  #prediction of each label) . I have to transpose one or the other (all_theta or X...which one? 
  #as i need to get a vector of probabilities for each label after multiplying them)
  
  #this will have the prediction prob in different col for a row from X. taking a max row wise will give
  #us the required p 
  pred_mat <- H(X %*% t(all_theta))
  
  p <- apply(pred_mat, 1, which.max)
  p  
}
