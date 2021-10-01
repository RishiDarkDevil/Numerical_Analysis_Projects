library(rgl)
library(plot3D)
#library(misc3d)
#library(akima)
library(tidyverse)
library(plotly)

# Number of partitions on each side
(n <- 50)

# Creating the boundaries of the square
(edge.1 <- rep(0, n+1))
(edge.2 <- seq(0, 1, length.out = n+1))
(edge.3 <- edge.1 + 1)
(edge.4 <- edge.2^3)

# List of all the boundary heights
c2 <- c(edge.3)
for (i in 2:n) {
  c2 <- c(c2, rev(edge.4)[i], rev(edge.2)[i])
}
c2 <- c(c2, edge.1)
c2
length(c2)
# |-----Edge.3------|
# |                 |
# |                 |
# |Edge.4           |Edge.2
# |                 |
# |                 |
# |-----Edge.1------|


(z <- c(edge.1, edge.4, edge.3, edge.2))
(x <- c(edge.2, edge.3, rev(edge.2), edge.1))
(y <- c(edge.1, edge.2, edge.3, edge.2))

bd.frame <- tibble(x = x, y = y, z = z)

# ---------------------- +ve X-axis
# |
# |    Z-axis is outwards
# |
# | -ve Y-axis

scatter3D(x = x, y = y, z = z)
points3d(x = x, y = y, z = z)

# Initializing some important stuffz
(gap <- 1/n)
(T.area <- 0.5 * gap * gap) # Area of the triangle
(nbdpts <- length(c2))
(pts <- (n+1)*(n+1)) # All the points
(nt <- n*n*2) # Number of triangles


# ---------------------- +ve X-axis
# |1 \ 2| 3\ 4| 5\ 6|7 \ 8| ---> My triangualtion is something like this  
# |9 \10| Z-axis is outwards
# |
# | -ve Y-axis

# A  F
# B 1 E --> 1 indicates a point and A,B,C,D,E,F refers to the triangles surrounding it in the resp. direction
#  C  D
# Let us code A-->1, B-->2, C-->3, D-->4, E-->5, F-->6

# Here I have changed the numbering pattern of the matrix a bit --> Interior points first are numbered to make partitioning easier
# 10  11  12  13  14
# 15  1   2   3   16
# 17  4   5   6   18
# 19  7   8   9   20
# 21  22  23  24  25


triangle.point.map <- matrix(rep(0, pts*nt), pts, nt) # Creating Points --> Triangle Numbers of which it is a part of
triangle.point.map
triangle.point.map.number <- matrix(rep(0, pts*6), nrow = pts, ncol = 6)
triangle.point.map.number

# Finding the edge points
edge.4.points <- c((n-1)^2 + 1)
edge.2.points <- c()
edge.1.points <- c((n-1)^2 + n + 2 + 2*(n-1))
edge.3.points <- c((n-1)^2 + 1)
for (i in 1:n){
  edge.4.points <- c(edge.4.points, (n-1)^2 + n + 2 + 2*(i-1))
  edge.2.points <- c(edge.2.points, (n-1)^2 + n + 1 + 2*(i-1))
  edge.1.points <- c(edge.1.points, (n-1)^2 + n + 2 + 2*(n-1) + i)
  edge.3.points <- c(edge.3.points, (n-1)^2 + 1 + i)
}
edge.2.points <- c(edge.2.points, pts)
edge.1.points
edge.2.points
edge.3.points
edge.4.points
bd.points <- c(edge.1.points, edge.2.points, edge.3.points, edge.4.points)
bd.points

# Filling triangle.point.map
for (i in 1:pts) {
  
  # Marking All the interior points
  if (!(i %in% bd.points)) {
    
    
    if (i %% (n-1) == 0) {
      triangle.point.map.number[i, 3] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2 # Marking C
      triangle.point.map.number[i, 4] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2 + 1 # Marking D
      triangle.point.map.number[i, 5] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2 + 1 + 1 # Marking E
      triangle.point.map.number[i, 2] = (2*i - 1) + (i%/%(n-1))*2 - 2 # Marking B
      triangle.point.map.number[i, 1] = (2*i - 1) + (i%/%(n-1))*2 - 2 + 1 # Marking A
      triangle.point.map.number[i, 6] = (2*i - 1) + (i%/%(n-1))*2 - 2 + 1 + 1 # Marking F
      
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2] = 3 # Marking C
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2 + 1] = 4 # Marking D
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 - 2 + 1 + 1] = 5 # Marking E
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2 - 2] = 2 # Marking B
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2 - 2 + 1] = 1 # Marking A
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2 - 2 + 1 + 1] = 6 # Marking F
    }
    else {
      triangle.point.map.number[i, 3] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 # Marking C
      triangle.point.map.number[i, 4] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 + 1 # Marking D
      triangle.point.map.number[i, 5] = (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 + 1 + 1 # Marking E
      triangle.point.map.number[i, 2] = (2*i - 1) + (i%/%(n-1))*2 # Marking B
      triangle.point.map.number[i, 1] = (2*i - 1) + (i%/%(n-1))*2 + 1 # Marking A
      triangle.point.map.number[i, 6] = (2*i - 1) + (i%/%(n-1))*2 + 1 + 1 # Marking F
      
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2] = 3 # Marking C
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 + 1] = 4 # Marking D
      triangle.point.map[i, (i%/%(n-1) + 1)*(n*2) + (i%%(n-1))*2 + 1 + 1] = 5 # Marking E
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2] = 2 # Marking B
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2 + 1] = 1 # Marking A
      triangle.point.map[i, (2*i - 1) + (i%/%(n-1))*2 + 1 + 1] = 6 # Marking F
    }
  }
  
  # Marking All the Edge.1 Points
  if (i %in% bd.points) {
    
    # Edge.1 Marking
    if (i %in% edge.1.points[-(n+1)]){
      triangle.point.map.number[i, 6] = 2*n*(n-1) + 2*(i%%(n+1)) - 1 # Marking F
      
      triangle.point.map[i, 2*n*(n-1) + 2*(i%%(n+1)) - 1] = 6 # Marking F
    }
    if (i %in% edge.1.points[-1]){
      triangle.point.map.number[i, 1] = 2*n*(n-1) + 2*(i%%(edge.1.points[1])) # Marking A
      triangle.point.map.number[i, 2] = 2*n*(n-1) + 2*(i%%(edge.1.points[1])) - 1 # Marking B
      
      triangle.point.map[i, 2*n*(n-1) + 2*(i%%(edge.1.points[1]))] = 1 # Marking A
      triangle.point.map[i, 2*n*(n-1) + 2*(i%%(edge.1.points[1])) - 1] = 2 # Marking B
    }
    
    
    # Edge.2 Marking
    # ---Done outside the end of this loop
    
    # Edge.3 Marking
    if (i %in% edge.3.points[-(n+1)]) {
      triangle.point.map.number[i, 5] = 2*(i%%(edge.3.points[1])) + 2 # Marking E
      triangle.point.map.number[i, 4] = 2*(i%%(edge.3.points[1])) + 2 - 1 # Marking D
      
      triangle.point.map[i, 2*(i%%(edge.3.points[1])) + 2] = 5 # Marking E
      triangle.point.map[i, 2*(i%%(edge.3.points[1])) + 2 - 1] = 4 # Marking D
    }
    if (i %in% edge.3.points[-c(1,n+1)]){
      triangle.point.map.number[i, 3] = 2*(i%%(edge.3.points[1])) + 2 - 1 - 1 # Marking C
      
      triangle.point.map[i, 2*(i%%(edge.3.points[1])) + 2 - 1 - 1] = 3 # Marking C
    }
    
    
    # Edge.4 Marking
    # ---Done outside the end of this loop
  }
}

# Edge.2 Marking
for (i in 1:n){
  triangle.point.map.number[edge.2.points[i], 3] = i*2*n # Marking C
  
  triangle.point.map[edge.2.points[i], i*2*n] = 3 # Marking C
}
for (i in 2:n){
  triangle.point.map.number[edge.2.points[i], 2] = (i-1)*2*n - 1 # Marking B
  triangle.point.map.number[edge.2.points[i], 1] = (i-1)*2*n # Marking A
  
  triangle.point.map[edge.2.points[i], (i-1)*2*n - 1] = 2 # Marking B
  triangle.point.map[edge.2.points[i], (i-1)*2*n] = 1 # Marking A
}

# Edge.4 Marking
for (i in 2:n) {
  triangle.point.map.number[edge.4.points[i], 4] = (i-1)*2*n # Marking D
  triangle.point.map.number[edge.4.points[i], 5] = (i-1)*2*n + 1 + 1 # Marking E
  triangle.point.map.number[edge.4.points[i], 6] = (i-1)*2*n + 1 - 2*n # Marking F
  
  triangle.point.map[edge.4.points[i], (i-1)*2*n + 1] = 4 # Marking D
  triangle.point.map[edge.4.points[i], (i-1)*2*n + 1 + 1] = 5 # Marking E
  triangle.point.map[edge.4.points[i], (i-1)*2*n + 1 - 2*n] = 6 # Marking F
}

triangle.point.map

triangle.point.map.number

# The slopes of A, B, C, D, E, F planes 
beta.ij <- c(0, 1/gap, 1/gap, 0, -1/gap, -1/gap)
gamma.ij <- c(-1/gap, 0, 1/gap, 1/gap, 0, -1/gap)

beta.slope <- function(k) {
  if (k == 0)
    return(beta.ij[1])
  else
    return(beta.ij[k])
} 
gamma.slope <- function(k) {
  if (k == 0)
    return(0)
  else
    return(gamma.ij[k])
} 

# This is our M Matrix
M <- matrix(rep(0, pts*(pts - length(c2))), nrow = (pts - length(c2)))

# Saving on Computation
triangle.point.map.beta <- matrix(rep(0, pts*nt), pts, nt)
triangle.point.map.gamma <- matrix(rep(0, pts*nt), pts, nt)

# Trying to improve computation speed
for (i in 1:pts){
  triangle.point.map.beta[i,triangle.point.map.number[i,]] <- unlist(lapply(triangle.point.map[i,triangle.point.map.number[i,]], beta.slope))
  triangle.point.map.gamma[i,triangle.point.map.number[i,]] <- unlist(lapply(triangle.point.map[i,triangle.point.map.number[i,]], gamma.slope))
}
triangle.point.map.beta
triangle.point.map.gamma

# Filling M -- !!!! WARNING !!! Takes lot of time!!!
for (i in 1:(pts - length(c2))) {
  for (j in 1:pts) {
    if (abs(i-j) <= 1 | abs(i-j) == n-1 | abs(i-j) == n | i %in% bd.points | j %in% bd.points) {
      M[i, j] = T.area*sum(triangle.point.map.beta[i,triangle.point.map.number[i,]]*triangle.point.map.beta[j,triangle.point.map.number[i,]] + triangle.point.map.gamma[i,triangle.point.map.number[i,]]*triangle.point.map.gamma[j,triangle.point.map.number[i,]])
    }
  }
}

M

# Trying to Solve
M11 <- M[1:(pts - length(c2)), 1:(pts - length(c2))]
M11
M12 <- M[1:(pts - length(c2)), (pts - length(c2) + 1):pts]
M12
RHS <- -M12 %*% c2
RHS
c1 <- solve(M11, RHS)
c1 <- Solve.Gauss.Seidel(M11, RHS, seq(0, 1, length.out = length(RHS)), 1000, 1e-5)
c1 <- Solve.Gauss.Jacobi(M11, RHS, seq(0, 1, length.out = length(RHS)), 1000, 1e-3)
c1
# Trying to Plot
(x1 <- seq(0+gap, 1-gap, length.out = n-1))
(y1 <- rev(seq(0+gap, 1-gap, length.out = n-1)))
soap.film <- cbind(expand.grid(x1, y1), c1)
colnames(soap.film) <- c('x', 'y', 'z')
initial <- cbind(x, y, z)
soap.film <- rbind(initial, soap.film)
soap.film
points3d(x = soap.film$x, y = soap.film$y, z = soap.film$z)

# Trying 3D Plot
(x1 <- seq(0+gap, 1-gap, length.out = n-1))
(y1 <- rev(seq(0+gap, 1-gap, length.out = n-1)))
soap.film <- cbind(expand.grid(x1, y1), c1)
colnames(soap.film) <- c('x', 'y', 'z')
soap.film
soap.film <- as_tibble(soap.film) %>%
  arrange(x,y)
soap.film

z.matrix <- matrix(rep(0, (n-1)^2), n-1)
k <- 1
for (i in 1:(n-1)) {
  for (j in 1:(n-1)) {
    z.matrix[i, j] <- soap.film$z[k]
    k = k+1
  }
}

z.matrix

a <- mesh(x1,y1)
u <- a$x ; v <- a$y
x <- u
y <- v
z <- soap.film$z
surf3D(x, y, z.matrix, colvar = z.matrix, colkey = TRUE, 
       box = TRUE, bty = "b", phi = 20, theta = 120)

col <- rainbow(10)[cut(z.matrix, breaks = 10)]
surface3d(x, y, z.matrix, color = col)
persp3d(x, y, z.matrix, color = col)

z <- z.matrix

plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~z.matrix)

int.points <- tibble(x = x, y = y)

fig <- plot_ly(bd.frame %>% mutate(x = 1-x, y = 1-y), x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines+markers',
               line = list(width = 6, color = ~c, colorscale = 'Viridis'),
               marker = list(size = 3.5, color = ~c, colorscale = 'Greens'))
fig %>% add_surface(data = int.points, x = ~x, y = ~y, z = ~z)


# Trying 3D Plot
(x1 <- seq(0, 1, length.out = n+1))
(y1 <- rev(seq(0, 1, length.out = n+1)))

q <- c(rev(edge.3))
for(i in 2:n){
  q <- c(q, rev(edge.4)[i],c1[((n-1)*(i-2)+1):((n-1)*(i-2)+1+(n-2))], rev(edge.2)[i])
}
(q <- c(q, rev(edge.1)))

soap.film <- cbind(expand.grid(x1, y1), q)
colnames(soap.film) <- c('x', 'y', 'z')
head(soap.film)
soap.film <- as_tibble(soap.film) %>%
  arrange(x,y)
soap.film

z.matrix <- matrix(rep(0, (n+1)^2), n+1)
k <- 1
for (i in 1:(n+1)) {
  for (j in 1:(n+1)) {
    z.matrix[i, j] <- soap.film$z[k]
    k = k+1
  }
}

z.matrix

a <- mesh(x1,y1)
u <- a$x ; v <- a$y
x <- u
y <- v
z <- soap.film$z
#surf3D(x, y, z.matrix, colvar = z.matrix, colkey = TRUE, 
#       box = TRUE, bty = "b", phi = 20, theta = 120)

#col <- rainbow(10)[cut(z.matrix, breaks = 10)]
#surface3d(x, y, z.matrix, color = col)
#persp3d(x, y, z.matrix, color = col)

z <- z.matrix
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~z.matrix)

all.points <- tibble(x = x, y = y)

fig <- plot_ly(bd.frame %>% mutate(x = 1-x, y = 1-y), x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines+markers',
               line = list(width = 15, color = '#DA16FF'),
               marker = list(size = 5, color = '#AA0DFE'))
fig %>% add_surface(data = all.points, x = ~x, y = ~y, z = ~z)


# Basis Functions
# Trying 3D Plot
(x1 <- seq(0, 1, length.out = n+1))
(y1 <- rev(seq(0, 1, length.out = n+1)))


soap.film <- cbind(expand.grid(x1, y1), rep(0, nrow(expand.grid(x1, y1))))
soap.film[n,3] <- 1
soap.film[(n-2)^2,3] <- 1
colnames(soap.film) <- c('x', 'y', 'z')
head(soap.film)
soap.film <- as_tibble(soap.film) %>%
  arrange(x,y)
soap.film

z.matrix <- matrix(rep(0, (n+1)^2), n+1)
k <- 1
for (i in 1:(n+1)) {
  for (j in 1:(n+1)) {
    z.matrix[i, j] <- soap.film$z[k]
    k = k+1
  }
}

z.matrix

a <- mesh(x1,y1)
u <- a$x ; v <- a$y
x <- u
y <- v
z <- soap.film$z
#surf3D(x, y, z.matrix, colvar = z.matrix, colkey = TRUE, 
#       box = TRUE, bty = "b", phi = 20, theta = 120)

#col <- rainbow(10)[cut(z.matrix, breaks = 10)]
#surface3d(x, y, z.matrix, color = col)
#persp3d(x, y, z.matrix, color = col)

z <- z.matrix
plot_ly() %>%
  add_surface(x = ~x, y = ~y, z = ~z.matrix)

all.points <- tibble(x = x, y = y)

fig <- plot_ly(bd.frame %>% mutate(x = 1-x, y = 1-y), x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines+markers',
               line = list(width = 15, color = '#DA16FF'),
               marker = list(size = 5, color = '#AA0DFE'))
fig %>% add_surface(data = all.points, x = ~x, y = ~y, z = ~z)