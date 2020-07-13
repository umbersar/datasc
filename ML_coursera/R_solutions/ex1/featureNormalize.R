
meanNormalizeUsingSD <- function(x){
    #if sd is 0, then return the vector as it is. It would happen for dummy feature
    # ifelse(test = sd(x)!=0, yes = (x- mean(x))/(sd(x)), no = x)
  if(sd(x)!=0) (x- mean(x))/(sd(x)) else x
}

featureNormalize <- function(X) {

  X_norm <- X
  mu <- rep(0,ncol(X))
  sigma <- rep(0,ncol(X))
  
  mu = apply(X,2,mean)
  sigma = apply(X,2,sd)
  
  # X_norm = apply(X, 2, meanNormalizeUsingSD)#no need to apply as we vectorize. see below 
  (X-mu)/(sigma)#this is vectorized code. Thus no need for apply
  
  list(X_norm = X_norm, mu = mu, sigma = sigma)
  # ------------------------------------------------------------
}
