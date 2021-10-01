# Finding Bases for Column Space, Row Space and Null Space of a Matrix

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

# Gaussian Elimination for Reduced Row Echelon Form
rrefmatrix <- function(m){
  count.rows <- nrow(m)
  count.cols <- ncol(m)
  rank.count <- 0
  piv.cols <- c()
  piv <- 1
  H <- diag(1, count.rows)
  
  for (row.curr in 1:count.rows) {
    if(piv <= count.cols){
      i <- row.curr
      while (m[i,piv] == 0 && i <= count.rows) {
        i <- i+1
        if(i>count.rows){
          i <- row.curr
          piv <- piv+1
          if (piv > count.cols){
            return(list('RREF' = m, 'Rank' = rank.count, 'Pivs' = piv.cols, 'H' = H))
          }
        }
      }
      rank.count <- rank.count + 1
      if (i != row.curr) {
        m <- swaprows(m, i, row.curr)
        rownames(m)[i] <- row.curr
        rownames(m)[row.curr] <- i
        H1 <- diag(1, count.rows)
        H1 <- swaprows(H1, i, row.curr)
        H <- H1 %*% H
      }
      p <- 1/m[row.curr,piv]
      m <- scalerow(m, row.curr, 1/m[row.curr,piv])
      H2 <- diag(1, count.rows)
      H2 <- scalerow(H2, row.curr, p)
      H <- H2 %*% H
      for (j in 1:count.rows)
        if (j != row.curr) {
          k = m[j,piv] / m[row.curr,piv]
          m <- replacerow(m, row.curr, j, -k)
          H3 <- diag(1, count.rows)
          H3 <- replacerow(H3, row.curr, j, -k)
          H <- H3 %*% H
        }
      piv.cols <- append(piv.cols, piv)
      piv <- piv+1
    }
  }
  return(list('RREF' = m, 'Rank' = rank.count, 'Pivs' = piv.cols, 'H' = H))
}

# Finding Bases of Column Space, Rowspace and Null Space of a Matrix
Find.Bases <- function(A){
  rrefA <- rrefmatrix(A)
  r <- rrefA$Rank
  
  # Finding the Column Space Basis
  print("Column Space Basis Vectors:")
  for (i in rrefA$Pivs){
    print(A[,i])
  }
  # Finding the Row Space Basis
  print("Row Space Basis Vectors:")
  for (i in 1:min(nrow(A), ncol((A)))){
    if(rrefA$RREF[i,i] != 0)
      print(rrefA$RREF[i,])
  }
  # Finding the Null Space Basis
  print("Null Space Basis Vectors:")
  if(r == ncol(A)){
    print("NULL")
  }
  else{
    H <- rrefmatrix(t(rrefA$RREF))$H
    for (i in (nrow(H)-(ncol(A)-r)+1):nrow(H))
      print(H[i,])
  }
}

# Testing
mat <- matrix(1:12, 4, 3)
print(mat)
rrefmatrix(mat)
Find.Bases(mat)

m <- matrix(c(1,1,1,2,2,2,3,3,3), 3)
print(m)
rrefmatrix(m)
Find.Bases(m)

m <- matrix(sample(100, 50, replace = T), 5, 10)
print(m)
rrefmatrix(m)
Find.Bases(m)

