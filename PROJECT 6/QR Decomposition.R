library(pracma)
library(Matrix)

# QR Decomposition

# Efficient Implementation

# Detecting if input matrix is full column rank
rrefmatrix <- function(m){
  count.rows <- nrow(m)
  count.cols <- ncol(m)
  rank.count <- 0
  piv.cols <- c()
  piv <- 1
  
  for (row.curr in 1:count.rows) {
    if(piv <= count.cols){
      i <- row.curr
      while (m[i,piv] == 0 && i <= count.rows) {
        i <- i+1
        if(i>count.rows){
          i <- row.curr
          piv <- piv+1
          if (piv > count.cols){
            return(list('RREF' = m, 'Rank' = rank.count, 'Pivs' = piv.cols))
          }
        }
      }
      rank.count <- rank.count + 1
      if (i != row.curr) {
        m <- swaprows(m, i, row.curr)
        rownames(m)[i] <- row.curr
        rownames(m)[row.curr] <- i
      }
      m <- scalerow(m, row.curr, 1/m[row.curr,piv])
      for (j in 1:count.rows)
        if (j != row.curr) {
          k = m[j,piv] / m[row.curr,piv]
          m <- replacerow(m, row.curr, j, -k)
        }
      piv.cols <- append(piv.cols, piv)
      piv <- piv+1
    }
  }
  return(list('RREF' = m, 'Rank' = rank.count, 'Pivs' = piv.cols))
}

# Making Unit vectors
unit <- function(v) {v/sqrt(sum(v*v))}

# Doing multiplication with Householder Matrix
hmult <- function(u, x) {x - 2*sum(u*x)*u}

# Creating shaver for shaving off the cols to concentrate norm on top 
shaver <- function(x){
  x[1] <- x[1] - sqrt(sum(x*x))
  unit(x)
}

# Computing R from the packed QR Decomposed Matrix
compute.R <- function(m,R.diag){
  for (i in 1:ncol(m)) {
    m[i:nrow(m),i] = numeric(nrow(m)-i+1)
    m[i,i] = R.diag[i]
  }
  return(list("R" = m))
}

# Computing Qx when the packed QR Decomposed Matrix is passed along with x vector
multiply.Q <- function(m,x){
  res.vect <- x
  ulist <- list()
  for (i in ncol(m):1) {
    u = m[i:nrow(m),i]
    v <- numeric(nrow(m))
    v[i:nrow(m)] = u
    name <- paste("u",i,":",sep = "")
    ulist[[name]] <- v
    res.vect <- hmult(v,res.vect)
  }
  return(list("Qx"=res.vect, "ulist"=ulist))
}

# Computing Q from the packed QR Decomposed Matrix - Not Needed Did extra
compute.Q <- function(m){
  e <- numeric(nrow(m))
  e[1] = 1
  A <- multiply.Q(m,e)
  Q <- A$Qx
  
  for (i in 2:nrow(m)){
    e <- numeric(nrow(m))
    e[i] = 1
    Q <- rbind(Q,multiply.Q(m,e)$Qx)
  }
  return(list("Q" = t(Q), "U" = rev(A$ulist)))
}

# Efficient QR decomposition
Eff.QR <- function(mat){
  R.diags = numeric(ncol(mat))
  rref <- rrefmatrix(mat)
  pivs <- rref$Pivs
  if(rref$Rank != ncol(mat))
    warning("Matrix is not Full Column Rank")
  for (i in 1:ncol(mat)) {
    if(i %in% pivs)
      u = shaver(mat[i:nrow(mat),i])
    else
      u = rep(0, nrow(mat)-i+1)
    for (j in i:ncol(mat)) {
      mat[i:nrow(mat),j] <- hmult(u, mat[i:nrow(mat),j])
    }
    if(i %in% pivs)
      R.diags[i] = mat[i,i]
    else
      R.diags[i] = 0
    mat[i:nrow(mat),i] = u
  }
  computeQ <- compute.Q(mat)
  return(list("Packed_Matrix" = mat, "R_Diag" = R.diags, "Q" = computeQ$Q, "R" = compute.R(mat, R.diags)$R, "U" = computeQ$U))
}

#Testing: Successful
mat <- matrix(c(1,1,1,7,-2,2,3,-2,2,8,3,2,5,1,5,4,4,0,3,4), nrow = 5, ncol = 4)
print(mat)
Eff.QR(mat)

 # Doing some extra work: Computing Q
compute.Q <- function(m){
  e <- numeric(nrow(m))
  e[1] = 1
  A <- multiply.Q(m,e)
  Q <- A$Qx
  
  for (i in 2:nrow(m)){
    e <- numeric(nrow(m))
    e[i] = 1
    Q <- rbind(Q,multiply.Q(m,e)$Qx)
  }
  return(list("Q" = t(Q), "U" = rev(A$ulist)))
}

compute.Q(m)

Q <- diag(1, nrow(m), nrow(m)) 

for (i in 1:ncol(m)) {
  u = matrix(m[i:nrow(m),i])
  H = diag(1, nrow = length(u), ncol = length(u)) - 2*u%*%t(u)
  actual.H <- diag(1,nrow(m),nrow(m))
  actual.H[i:nrow(actual.H),i:ncol(actual.H)] <-  H
  Q <- Q%*%actual.H
}

print(Q)

# Doing some extra work: Computing R
compute.R <- function(m){
  for (i in 1:ncol(m)) {
    m[i:nrow(m),i] = numeric(nrow(m)-i+1)
    m[i,i] = R.diag[i]
  }
  return(list("R" = m))
}

compute.R(m)

Q%*%R # For checking if I get back my original matrix
Q %*% t(Q) # For checking Orthogonality

outMatrix <- matrix(data = 1:12,4,3)
print(outMatrix)
Eff.QR(outMatrix)
multiply.Q(Eff.QR(outMatrix)$Packed_Matrix, c(4,1,5,0))$Qx

mat <- matrix(c(1,1,1,2,2,2,3,3,3), nrow = 3, ncol = 3)
print(mat)
Eff.QR(mat)


# Computing Qx where x is a vector in a efficient way
x <- c(4,1,5,0)

multiply.Q <- function(m,x){
  res.vect <- x
  ulist <- list()
  for (i in ncol(m):1) {
    u = m[i:nrow(m),i]
    v <- numeric(nrow(m))
    v[i:nrow(m)] = u
    name <- paste("u",i,":",sep = "")
    ulist[[name]] <- v
    res.vect <- hmult(v,res.vect)
  }
  return(list("Qx"=res.vect, "ulist"=ulist))
}

Q %*% x

householder(outMatrix) # Checking with Pracma package


# Downloaded from Internet rpubs.com

qr.householder <- function(A) {
  require(Matrix)
  
  R <- as.matrix(A) # Set R to the input matrix A
  
  n <- ncol(A)
  m <- nrow(A)
  H <- list() # Initialize a list to store the computed H matrices to calculate Q later
  
  if (m > n) {
    c <- n
  } else {
    c <- m
  }
  
  for (k in 1:c) {
    x <- R[k:m,k] # Equivalent to a_1
    e <- as.matrix(c(1, rep(0, length(x)-1)))
    vk <- sign(x[1]) * sqrt(sum(x^2)) * e + x
    
    # Compute the H matrix
    hk <- diag(length(x)) - 2 * as.vector(vk %*% t(vk)) / (t(vk) %*% vk)
    if (k > 1) {
      hk <- bdiag(diag(k-1), hk)
    }
    
    # Store the H matrix to find Q at the end of iteration
    H[[k]] <- hk
    
    R <- hk %*% R
  }
  
  Q <- Reduce("%*%", H) # Calculate Q matrix by multiplying all H matrices
  res <- list('Q'=Q,'R'=R)
  return(res)
}

qr.householder(outMatrix)




(A <- matrix(c(2,8,0,1,0,0,0,8,0), 3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:3, sep = '')
rrefmatrix(A)
