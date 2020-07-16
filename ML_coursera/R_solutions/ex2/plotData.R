library(ggplot2)
plotData <-
  function (X, y, axLables = c("Exam 1 score","Exam 2 score"), legLabels =
              c('Admitted', 'Not admitted')) {
    #PLOTDATA Plots the data points X and y into a new device
    #   PLOTDATA(x,y) plots the data points with + for the positive examples
    #   and o for the negative examples. X is assumed to be a Mx2 matrix.
    
    # ----------------------- YOUR CODE HERE -----------------------
    # Instructions: Plot the positive and negative examples on a
    #               2D plot, using the option pch=3 for the positive (plus)
    #               examples and pch=21 for the negative (circle) examples.
    #
    
    df <- data.frame(X,y)
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], color=y)) + geom_point(shape=1)
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], color=factor(y))) + geom_point(shape=1)
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], color=factor(y))) + geom_point()
    ggplot(df, aes(x=X[,1], y=X[,2])) +geom_point(aes(color = factor(y)))
    
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], shape=y)) + geom_point() + scale_shape_identity()
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], shape=factor(y))) + geom_point() + scale_shape_identity()
    ggplot2::ggplot(df, aes(x=X[,1], y=X[,2], shape=factor(y))) + geom_point() 
    # ----------------------------------------------------
}