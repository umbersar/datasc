
meanNormalizeUsingSD <- function(x, mu,sigma){
    #if sd is 0, then return the vector as it is. It would happen for dummy feature
    ifelse(test = sd(x)!=0, yes = (x- mean(x))/(sd(x)), no = x)
}

featureNormalize <- function(X) {

  X_norm <- X
  mu <- rep(0,ncol(X))
  sigma <- rep(0,ncol(X))
  
  mu = apply(X,2,mean)
  sigma = apply(X,2,sd)
  
  X_norm = apply(X, 2, meanNormalizeUsingSD)
  
  list(X_norm = X_norm, mu = mu, sigma = sigma)
  # ------------------------------------------------------------
}
