# GRAVITY-WELL
library(plot3D)
library(rgl)
library(plotly)
library(processx)

g = 9.8

u = function(x, y) sqrt(x^2+y^2)
f = function(x, y) sqrt(u(x,y) - 1)
f.dash = function(x, y) 1/(2*sqrt(u(x,y) - 1))
R. <- function(x, y, vx, vy) ( ( f.dash(x, y) * (vx^2 + vy^2 - ((x*vx+y*vy)/u(x,y))^2 ) / u(x,y) ) + (((x*vx+y*vy)/u(x,y))^2)*(-1/(4*(u(x,y)-1)^1.5)) + g)/((u(x,y))*(f.dash(x, y) + (1/f.dash(x, y))))
a.x <- function(x, y, vx, vy) -x*R.(x, y, vx, vy)
a.y <- function(x, y, vx, vy) -y*R.(x, y, vx, vy)

ts <- function(xyz,sleep=0.3){
  plot3d(xyz,type="n")
  n = nrow(xyz)
  p = points3d(xyz[1,])
  l = lines3d(xyz[1,])
  for(i in 2:n){
    Sys.sleep(sleep)
    rgl.pop("shapes",p)
    rgl.pop("shapes",l)
    p=points3d(xyz[i,])
    l=lines3d(xyz[1:i,])
  }
}

gravity.well = function(t0, x0, y0, vx0, vy0, n, dt){
  x = rep(0, n)
  y = rep(0, n)
  z = rep(0, n)
  t = rep(0, n)
  vx = rep(0,n)
  vy = rep(0,n)
  # ax = rep(0,n)
  # ay = rep(0,n)
  
  x[1] = x0
  y[1] = y0
  z[1] = f(x[1], y[1])
  t[1] = t0
  vx[1] = vx0
  vy[1] = vy0
  # ax[1] = 0
  # ay[1] = 0
  h <- dt
  
  for (i in 2:n) {
    t[i] = t[i-1] + dt
    k1 <- h*vx[i-1]
    j1 <- h*vy[i-1]
    l1 <- h*a.x(x[i-1], y[i-1], vx[i-1], vy[i-1])
    m1 <- h*a.y(x[i-1], y[i-1], vx[i-1], vy[i-1])
    k2 <- h*(vx[i-1] + l1/2)
    j2 <- h*(vy[i-1] + m1/2)
    l2 <- h*a.x(x[i-1] + k1/2, y[i-1] + j1/2, vx[i-1] + l1/2, vy[i-1] + m1/2)
    m2 <- h*a.y(x[i-1] + k1/2, y[i-1] + j1/2, vx[i-1] + l1/2, vy[i-1] + m1/2)
    k3 <- h*(vx[i-1] + l2/2)
    j3 <- h*(vy[i-1] + m2/2)
    l3 <- h*a.x(x[i-1] + k2/2, y[i-1] + j2/2, vx[i-1] + l2/2, vy[i-1] + m2/2)
    m3 <- h*a.y(x[i-1] + k2/2, y[i-1] + j2/2, vx[i-1] + l2/2, vy[i-1] + m2/2)
    k4 <- h*(vx[i-1] + l3)
    j4 <- h*(vy[i-1] + m3)
    l4 <- h*a.x(x[i-1] + k3, y[i-1] + j3, vx[i-1] + l3, vy[i-1] + m3)
    m4 <- h*a.y(x[i-1] + k3, y[i-1] + j3, vx[i-1] + l3, vy[i-1] + m3)
    
    vx[i] = vx[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6 
    vy[i] = vy[i-1] + (m1 + 2*m2 + 2*m3 + m4)/6
    x[i] = x[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6
    y[i] = y[i-1] + (j1 + 2*j2 + 2*j3 + j4)/6
    z[i] = f(x[i], y[i])
    
    # Sys.sleep(0.01)
  }
  
  X = seq(-18, 18, 0.1)
  Y = seq(-18, 18, 0.1)
  z.matrix <- matrix(rep(0, length(X)*length(Y)), length(X))
  
  for (i in 1:length(X)) {
    for (j in 1:length(Y)) {
      z.matrix[i, j] <- f(X[i], Y[j])
    }
  }
  
  z.matrix[is.na(z.matrix)] <- 0
  
  for (i in 1:length(X)) {
    for (j in 1:length(Y)) {
      if(z.matrix[i, j] < 2.5 | z.matrix[i, j] == 0)
        z.matrix[i, j] = NA
    }
  }
  Z <- z.matrix
  
  # Plotting
  p <- plot_ly(source = "plot2") %>%
    add_surface(x = ~X, y = ~Y, z = ~Z) %>%
    add_markers(x = ~x, y = ~y, z = ~z, opacity = 0.8, color = c("red"), size = 0.5)
  return(p)
  
}

gravity.well(0, 10, 0, 0, 5, 50000, 0.01)
gravity.well(0, 5, 1, 0, 5, 50000, 0.1)
X = seq(-10, 10, 0.1)
Y = seq(-10, 10, 0.1)
Z = rep(0,200)
for (i in 1:201) {
  Z[i] = f(X[i], Y[i])
}

scatter3D(x = X, y = Y, z = Z, xlim = c(-10, 10), ylim = c(-10,10), zlim = c(0,10))

# -- Animation
gravity.well = function(t0, x0, y0, vx0, vy0, n, dt){
  x = rep(0, n)
  y = rep(0, n)
  z = rep(0, n)
  t = rep(0, n)
  vx = rep(0,n)
  vy = rep(0,n)
  # ax = rep(0,n)
  # ay = rep(0,n)
  
  x[1] = x0
  y[1] = y0
  z[1] = f(x[1], y[1])
  t[1] = t0
  vx[1] = vx0
  vy[1] = vy0
  # ax[1] = 0
  # ay[1] = 0
  h <- dt
  
  for (i in 2:n) {
    t[i] = t[i-1] + dt
    k1 <- h*vx[i-1]
    j1 <- h*vy[i-1]
    l1 <- h*a.x(x[i-1], y[i-1], vx[i-1], vy[i-1])
    m1 <- h*a.y(x[i-1], y[i-1], vx[i-1], vy[i-1])
    k2 <- h*(vx[i-1] + l1/2)
    j2 <- h*(vy[i-1] + m1/2)
    l2 <- h*a.x(x[i-1] + k1/2, y[i-1] + j1/2, vx[i-1] + l1/2, vy[i-1] + m1/2)
    m2 <- h*a.y(x[i-1] + k1/2, y[i-1] + j1/2, vx[i-1] + l1/2, vy[i-1] + m1/2)
    k3 <- h*(vx[i-1] + l2/2)
    j3 <- h*(vy[i-1] + m2/2)
    l3 <- h*a.x(x[i-1] + k2/2, y[i-1] + j2/2, vx[i-1] + l2/2, vy[i-1] + m2/2)
    m3 <- h*a.y(x[i-1] + k2/2, y[i-1] + j2/2, vx[i-1] + l2/2, vy[i-1] + m2/2)
    k4 <- h*(vx[i-1] + l3)
    j4 <- h*(vy[i-1] + m3)
    l4 <- h*a.x(x[i-1] + k3, y[i-1] + j3, vx[i-1] + l3, vy[i-1] + m3)
    m4 <- h*a.y(x[i-1] + k3, y[i-1] + j3, vx[i-1] + l3, vy[i-1] + m3)
    
    vx[i] = vx[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6 
    vy[i] = vy[i-1] + (m1 + 2*m2 + 2*m3 + m4)/6
    x[i] = x[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6
    y[i] = y[i-1] + (j1 + 2*j2 + 2*j3 + j4)/6
    z[i] = f(x[i], y[i])
    
    # Sys.sleep(0.01)
  }
  
  X = seq(-18, 18, 0.1)
  Y = seq(-18, 18, 0.1)
  z.matrix <- matrix(rep(0, length(X)*length(Y)), length(X))
  
  for (i in 1:length(X)) {
    for (j in 1:length(Y)) {
      z.matrix[i, j] <- f(X[i], Y[j])
    }
  }
  
  z.matrix[is.na(z.matrix)] <- 0
  
  for (i in 1:length(X)) {
    for (j in 1:length(Y)) {
      if(z.matrix[i, j] < 2.5 | z.matrix[i, j] == 0)
        z.matrix[i, j] = NA
    }
  }
  Z <- z.matrix
  
  # Plotting
  for (i in seq(1, n, by = 50)) {
    p <- plot_ly(source = "plot2") %>%
      add_surface(x = ~X, y = ~Y, z = ~Z) %>%
      add_markers(x = ~x[1:i], y = ~y[1:i], z = ~z[1:i], opacity = 0.8, color = c("red"), size = 2)
    orca(p, paste0("p",i,".png"))
  }
}

gravity.well(0, 10, 0, 0, 5, 5000, 0.01)

png_files <- list.files(getwd(), pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 0.1)
