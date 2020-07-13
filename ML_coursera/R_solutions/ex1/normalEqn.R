normalEqn <- function(X, y) {
  source("pinv.R")
  theta <- rep(0,length(y))
  
  theta <- pinv(t(X)%*%X) %*% t(X) %*% y
}
