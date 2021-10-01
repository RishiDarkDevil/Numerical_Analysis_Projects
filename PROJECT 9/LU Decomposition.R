# Computes the L matrix 
L <- function(A, j) {
  n <- nrow(A)
  tot <- rep(0, n-j+1)
  if (j > 1){
    for (k in 1:(j-1))
      tot <- tot + A[j:n, k] * A[k, j]
  }
  A[(j:n), j] <- A[(j:n), j] - tot
  return(A)
}

# Computes the U matrix
U <- function(A, i) {
  n <- nrow(A)
  tot <- rep(0, n-i)
  if (i > 1){
    for (k in 1:(i-1))
      tot <- tot + A[k, (i+1):n] * A[i, k]
  }
  if(A[i, i] == 0)
    stop("Singular Matrix Detected!")
  A[i, (i+1):n] <- (A[i, (i+1):n] - tot) / A[i, i]
  return(A)
}

# Test Successful
A <- matrix(1:9, 3, 3); n <- 3
A
(A <- L(A, 1))
(A <- U(A, 1))
(A <- L(A, 2))
(A <- U(A, 2))
(A <- L(A, 3))

# Computes the LU Decomposition
Compute.LU <- function(A) {
  if(nrow(A) != ncol(A))
    stop("Matrix is not square!")
  
  n <- nrow(A)
  
  for(i in 1:n) {
    A <- L(A, i)
    if(i < n)
      A <- U(A, i)
  }
  
  L <- matrix(rep(0, n^2), n, n)
  U <- diag(n)
  
  for(j in 1:n)
    L[j:n,j] <- A[j:n, j]
  
  for(i in 1:(n-1))
    U[i, (i+1):n] <- A[i, (i+1):n]
  
  return(list("A" = A, "L"=L, "U"=U))
}

# Testing: Successful
(A <- matrix(c(1,2,3,5,6,7,9,11,12), 3, 3))

(decomp <- Compute.LU(A))
decomp$L %*% decomp$U

(A <- matrix(c(1,1,1,2,2,2,3,3,3), 3, 3))

(decomp <- Compute.LU(A))
decomp$L %*% decomp$U

# Solve Linear System
Solve.LU <- function(A, b){
  LU <- Compute.LU(A)
  n <- nrow(A)
  
  # Forward Substitution
  y <- rep(0, n)
  if (LU$A[1,1] == 0) 
    stop("No Solution!")
  y[1] <- b[1] / LU$A[1,1] 
  for (i in 2:n) {
    if (LU$A[i,i] == 0)
      stop("No Solution!")
    y[i] <- (b[i] - sum(LU$A[i, 1:(i-1)]*y[1:(i-1)])) / LU$A[i,i] 
  }
  
  # Backward Substitution
  x <- rep(0, n)
  x[n] <- y[n]
  for (i in (n-1):1) {
    x[i] <- y[i] - sum(LU$A[i, (i+1):n]*x[(i+1):n])
  }
  
  return(list("LU"=LU, "x"=x))
}

# Testing: Successful
(A <- matrix(c(2,4,-6,4,3,7,-10,6,0,2,0,4,0,0,1,5), 4, 4))

Compute.LU(A)
Compute.LU(A)$L %*% Compute.LU(A)$U

(sol <- Solve.LU(A, c(1,2,1,0)))
A %*% sol$x

