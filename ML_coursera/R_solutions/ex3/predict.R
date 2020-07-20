predict <- function(Theta1, Theta2, X) {
  #PREDICT Predict the label of an input given a trained neural network
  #   p <- PREDICT(Theta1, Theta2, X) outputs the predicted label of X given the
  #   trained weights of a neural network (Theta1, Theta2)
  
  #If there is only one example in the dataset which we have to predict, then transpose 
  #the feature vector to transform it into a feature matrix type shape with 1 row and n cols
  if (is.vector(X))
    X <- t(X)
  
  m <- nrow(X)#no. of rows in dataset for which to predict(features as columns and rows as examples)
  # add the dummy feature to the matrix of features(X...unlabeled data)
  X <- cbind(rep(1,m), X)
  
  num_labels <- nrow(Theta2)
  
  # vector containing labels between 1 to num_labels. return this
  p <- rep(0,m)
  
  z2 = X %*% t(Theta1)# we are applying the theta for all the different oneVsall classifiers to 
  #each example and the resulting matrix is one which have columns=no. of labels and rows=no of 
  #examples for which we have to predict
  a2 <- cbind(rep(1,nrow(z2)), H(z2))#add dummy feature for calculating next layer nodes/features
  
  z3 <- a2 %*% t(Theta2)
  #this will have the prediction probabilities in different cols for a row from a2. taking a max 
  #row wise will give us the required p 
  a3 <- H(z3)

  p <- apply(a3, 1, which.max)
  p
}
