#this is our hypothesis
H <- function(x) {
  1/(1 + exp(-x))
}

#this is just a wrapper to satisfy grader
sigmoid <- function(x) {
  H(x)
}