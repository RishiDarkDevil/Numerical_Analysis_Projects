# LEAST SQUARE SOLUTION to Ax = b

# Efficient Implementation

# Detecting if input matrix is full column rank

# Type-I Matrix
swaprows <- function(m, row1, row2){
  row.temp <- m[row1,]
  m[row1,] <- m[row2,]
  m[row2,] <- row.temp
  return(m)
}

# Type-II Matrix
scalerow <- function(m, row, k){
  m[row,] <- m[row,]*k
  return(m)
}

# Type-III Matrix
replacerow <- function(m, row1, row2, k){
  m[row2,] <- m[row2,] + m[row1,]*k
  return(m)
}

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

# Finding inverse using Gaussian Elimination
invmatrix <- function(m){
  count.rows <- nrow(m)
  count.cols <- ncol(m)
  piv <- 1
  inv <<- diag(1,count.rows)
  
  for (row.curr in 1:count.rows) {
    if(piv <= count.cols){
      i <- row.curr
      while (m[i,piv] == 0 && i <= count.rows) {
        i <- i+1
        if(i>count.rows){
          i <- row.curr
          piv <- piv+1
          if (piv > count.cols)
            return(m)
        }
      }
      if (i != row.curr) {
        m <- swaprows(m, i, row.curr)
        inv <<- swaprows(inv, i, row.curr)
        
      }
      p <- 1/m[row.curr,piv]
      m <- scalerow(m, row.curr, 1/m[row.curr,piv])
      inv <<- scalerow(inv, row.curr, p)
      
      for (j in 1:count.rows)
        if (j != row.curr) {
          k = m[j,piv] / m[row.curr,piv]
          m <- replacerow(m, row.curr, j, -k)
          inv <<- replacerow(inv, row.curr, j, -k)
        }
      piv <- piv+1
    }
  }
  return(list("Inverse" = inv))
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
multiply.Q <- function(m,x,t=1){ # t handles if we want to multiply x with Q? t = 1 implies Qx and t = 0 implies Q'x
  res.vect <- x
  ulist <- list()
  if(t == 1)
    k = ncol(m):1
  else
    k = 1:ncol(m)
  for (i in k) {
    u = m[i:nrow(m),i]
    v <- numeric(nrow(m))
    v[i:nrow(m)] = u
    name <- paste("u",i,":",sep = "")
    ulist[[name]] <- u
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
  if(rref$Rank != ncol(mat)){
    warning("Matrix is not Full Column Rank")
    stop("Can't Work with Arrays that is not Full Column Rank")
  }
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


# Finding Least Squares Solution
LS.solve <- function(A, b){
  QR <- Eff.QR(A)
  R1 <- QR$R[1:ncol(QR$R),]
  c1 <- multiply.Q(QR$Packed_Matrix, b, t=0)$Qx[1:ncol(QR$R)]
  R1.inv <- invmatrix(R1)$Inverse
  return(list("x" = R1.inv %*% c1))
}

# Finding Least Squares Solution - Better Version of the above function
LS.solve <- function(A, b){
  QR <- Eff.QR(A)
  R1 <- QR$R[1:ncol(QR$R),]
  c1 <- multiply.Q(QR$Packed_Matrix, b, t=0)$Qx[1:ncol(QR$R)]

  # Backward Substitution
  n <- nrow(R1)
  x <- rep(0, n)
  x[n] <- c1[n] / R1[n,n]
  for (i in (n-1):1) {
    x[i] <- (c1[i] - sum(R1[i, (i+1):n]*x[(i+1):n])) / R1[i,i]
  }
  
  return(list("x" = x))
}

# Testing
mat <- matrix(c(1,1,1,7,-2,2,3,-2,2,8,3,2,5,1,5,4,4,0,3,4), nrow = 5, ncol = 4)
print(mat)
Eff.QR(mat)
LS.solve(mat, c(1,2,3,4,5))
