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

# Gaussian Elimination to Row-Echelon Form
refmatrix <- function(m){
  count.rows <- nrow(m)
  count.cols <- ncol(m)
  piv <- 1
  
  for (row.curr in 1:count.rows) {
    if(piv <= count.cols){
      i <- row.curr
      while (m[i, piv] == 0 && i <= count.rows) {
        i <- i+1
        if(i > count.rows){
          i <- row.curr
          piv <- piv+1
          if(piv > count.cols)
            return(m)
        }
      }
      if(i != row.curr){
          m <- swaprows(m, i, row.curr)
          rownames(m)[i] <- row.curr
          rownames(m)[row.curr] <- i
      }
      for (j in row.curr:count.rows) 
        if (j != row.curr) {
          k <- m[j, piv] / m[row.curr, piv]
          m <- replacerow(m, row.curr, j, -k)
        }
      piv <- piv+1
    }
  }
  return(m)
}

(A <- matrix(c(2,8,5,4,4,5,0,8,4), 3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:3, sep = '')
refmatrix(A)

# Gaussian Elimination for Reduced Row Echelon Form
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

(A <- matrix(c(2,8,0,1,0,0,0,8,0), 3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:3, sep = '')
rrefmatrix(A)

(A <- matrix(c(2,8,0,1,0,0,0,8,0), 3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:3, sep = '')
rrefmatrix(A)

(A <- matrix(c(0,4,-1,3,1,-8,-5,1,-1, -1, 6, -4), nrow =3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:4, sep = '')
rrefmatrix(A)


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
}

(A <- matrix(c(2,0,0,-1,5,0,0,0,4), 3))
row.names(A) <- 1:3
colnames(A) <- paste("x", 1:3, sep = '')
invmatrix(A)
print(inv)
inv %*% A
diag(1,4)
scalerow(diag(1,4), 1, 0.5)







