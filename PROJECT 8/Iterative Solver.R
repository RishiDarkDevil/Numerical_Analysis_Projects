# Diagonally-Dominant Checker
#is.Diagonally.Dominant <- function(A){
#  for (i in 1:nrow(A)) {
#    if(abs(A[i,i]) <= (sum(abs(A[i,])) - abs(A[i,i])))
#      return(FALSE)
#  }
#  return(TRUE)
#}

# Checking Positive-Definite
#is.Positive.Definite <- function(A){
#  if(min(t(A) == A) == 0)
#    return(FALSE)
#  if(max(diag(Compute.LU(A)$L) < 0) == 1)
#    return(FALSE)
#  return(TRUE)
#}

# Gauss-Jacobi Method
Solve.Gauss.Jacobi <- function(A, b, guess, n, tol){
#  if (is.Diagonally.Dominant(A) == FALSE & is.Positive.Definite(A) == FALSE) {
#    stop("The System is neither Diagonally Dominant nor Positive Definite!")
#  }
  x <- guess
  dd <- diag(A)
  diag(A) <- rep(0, nrow(A)) 
  for (i in 1:n) {
    x_old <- x
    x <- (b - (A %*% x))/dd
    #print(t(x))
    if(sqrt(sum((x-x_old)^2)) < tol)
      return(x)
  }
  warning("Maximum Iteration Exceeded!")
  return(x)
}

# Testing: Successful
A <- matrix(c(20,1,1,3,-4,-4,-4,1,10), ncol = 3, nrow = 3)
b <- c(19, -3, 7)
guess <- c(0, 0, 0)
Solve.Gauss.Jacobi(A, b, guess, 10, 1e-5)

# Gauss-Seidel Method
Solve.Gauss.Seidel <- function(A, b, guess, n, tol){
#  if (is.Diagonally.Dominant(A) == FALSE & is.Positive.Definite(A) == FALSE) {
#    stop("The System is neither Diagonally Dominant nor Positive Definite!")
#  }
  x <- guess
  dd <- diag(A)
  diag(A) <- rep(0, nrow(A)) 
  for (i in 1:n) {
    x_old <- x
    for (k in 1:length(b))
      x[k] <- (b[k] - sum(A[k,]*x))/dd[k]
    #print(t(x))
    if(sqrt(sum((x-x_old)^2)) < tol)
      return(t(t(x)))
  }
  warning("Maximum Iteration Exceeded!")
  return(t(t(x)))
}

# Testing: Successful
A <- matrix(c(20,1,1,3,-4,-4,-4,1,10), ncol = 3, nrow = 3)
b <- c(19, -3, 7)
guess <- c(0, 0, 0)
Solve.Gauss.Seidel(A, b, guess, 10, 1e-5)
