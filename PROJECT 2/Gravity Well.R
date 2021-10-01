# GRAVITY-WELL
library(plot3D)
library(rgl)
library(plotly)
library(tidyverse)
library(gifski)

g = 9.8

u = function(x, y) sqrt(x^2+y^2)
f = function(x, y) sqrt(u(x,y) - 1)
f.dash = function(x, y) 1/(2*sqrt(u(x,y) - 1))

gravity.well = function(t0, x0, y0, vx0, vy0, n, dt){
  x = rep(0, n)
  y = rep(0, n)
  z = rep(0, n)
  t = rep(0, n)
  vx = rep(0,n)
  vy = rep(0,n)
  ax = rep(0,n)
  ay = rep(0,n)
  
  x[1] = x0
  y[1] = y0
  z[1] = f(x[1], y[1])
  t[1] = t0
  vx[1] = vx0
  vy[1] = vy0
  ax[1] = 0
  ay[1] = 0
  
  
  for (i in 2:n) {
    t[i] = t[i-1] + dt
    R = ( ( f.dash(x[i-1], y[i-1]) * (vx[i-1]^2 + vy[i-1]^2 - ((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2 ) / u(x[i-1],y[i-1]) ) + (((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2)*(-1/(4*(u(x[i-1],y[i-1])-1)^1.5)) + g)/((u(x[i-1],y[i-1]))*(f.dash(x[i-1], y[i-1]) + (1/f.dash(x[i-1], y[i-1])))) 
    ax[i] = -x[i-1]*R
    ay[i] = -y[i-1]*R
    vx[i] = vx[i-1] + ax[i-1]*dt 
    vy[i] = vy[i-1] + ay[i-1]*dt
    x[i] = x[i-1] + vx[i-1]*dt + ax[i-1]*(dt^2)/2
    y[i] = y[i-1] + vy[i-1]*dt + ay[i-1]*(dt^2)/2
    z[i] = f(x[i], y[i])
    
    # Sys.sleep(0.01)
  }
  
  #lines3D(x = x, y = y, z = z)
  #plot3d(x = x, y = y, z = z)
  
  X = seq(-23, 23, 0.1)
  Y = seq(-23, 23, 0.1)
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
  
  #col <- rainbow(10)[cut(z.matrix, breaks = 10)]
  #surface3d(X, Y, z.matrix, color = col)
  #plot3d(x = c(10), y = c(0), z = c(f(10, 0)), size = 15)
  Z <- z.matrix
  
  plot_ly() %>%
    add_surface(x = ~X, y = ~Y, z = ~Z) %>%
    add_markers(x = ~x, y = ~y, z = ~z, opacity = 0.8, color = c("red"), size = 0.5)
  
}

gravity.well(0, 10, 0, 0, 5, 50000, 0.01)
#X = seq(-10, 10, 0.1)
#Y = seq(-10, 10, 0.1)
#Z = rep(0,200)
#for (i in 1:201) {
#  Z[i] = f(X[i], Y[i])
#}

#scatter3D(x = X, y = Y, z = Z, xlim = c(-10, 10), ylim = c(-10,10), zlim = c(0,10))

X = seq(-23, 23, 0.1)
Y = seq(-23, 23, 0.1)
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

#col <- rainbow(10)[cut(z.matrix, breaks = 10)]
#surface3d(X, Y, z.matrix, color = col)
#plot3d(x = c(10), y = c(0), z = c(f(10, 0)), size = 15)
Z <- z.matrix
plot_ly() %>%
  add_surface(x = ~X, y = ~Y, z = ~Z) %>%
  add_markers(x = ~c(10), y = ~c(0), z = ~c(f(10, 0)))
surf3D(X,Y,z.matrix)

# ---- RGL ANIMATION

ts <- function(xyz,sleep=0.3){
  X = seq(-23, 23, 0.1)
  Y = seq(-23, 23, 0.1)
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
  col <- rainbow(10)[cut(z.matrix, breaks = 10)]
  
  plot3d(xyz,type="n", size = 10)
  surface3d(X, Y, z.matrix, color = col)
  n = nrow(xyz)
  p = points3d(xyz[1,], cex = 10)
  l = lines3d(xyz[1,], size = 10)
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
  ax = rep(0,n)
  ay = rep(0,n)
  
  x[1] = x0
  y[1] = y0
  z[1] = f(x[1], y[1])
  t[1] = t0
  vx[1] = vx0
  vy[1] = vy0
  ax[1] = 0
  ay[1] = 0
  
  
  for (i in 2:n) {
    t[i] = t[i-1] + dt
    R = ( ( f.dash(x[i-1], y[i-1]) * (vx[i-1]^2 + vy[i-1]^2 - ((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2 ) / u(x[i-1],y[i-1]) ) + (((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2)*(-1/(4*(u(x[i-1],y[i-1])-1)^1.5)) + g)/((u(x[i-1],y[i-1]))*(f.dash(x[i-1], y[i-1]) + (1/f.dash(x[i-1], y[i-1])))) 
    ax[i] = -x[i-1]*R
    ay[i] = -y[i-1]*R
    vx[i] = vx[i-1] + ax[i-1]*dt 
    vy[i] = vy[i-1] + ay[i-1]*dt
    x[i] = x[i-1] + vx[i-1]*dt + ax[i-1]*(dt^2)/2
    y[i] = y[i-1] + vy[i-1]*dt + ay[i-1]*(dt^2)/2
    z[i] = f(x[i], y[i])
    
    # Sys.sleep(0.01)
  }
  
  #lines3D(x = x, y = y, z = z)
  #plot3d(x = x, y = y, z = z)
  
  
  
  #col <- rainbow(10)[cut(z.matrix, breaks = 10)]
  #surface3d(X, Y, z.matrix, color = col)
  #plot3d(x = c(10), y = c(0), z = c(f(10, 0)), size = 15)
  xyz <- data.frame(x,y,z)
  ts(xyz, 1e-15)
  
  #plot_point <- tibble(t = t, x = x, y = y, z = z)
  
  #plot_ly(data = plot_point, x = ~x, y = ~y, z = ~z, opacity = 0.8, color = c("red"), size = 2, frame = ~t, mode = "markers", type = "scatter3d") %>%
    #layout(xaxis = list(range = c(-23, 23)), yaxis = list(range = c(-23, 23)))
    #add_surface(x = ~X, y = ~Y, z = ~Z)# %>%
    #add_markers(data = plot_point, x = ~x, y = ~y, z = ~z, opacity = 0.8, color = c("red"), size = 1, frame = ~t)
  
}

# --- ANIMATION
gravity.well = function(t0, x0, y0, vx0, vy0, n, dt){
  x = rep(0, n)
  y = rep(0, n)
  z = rep(0, n)
  t = rep(0, n)
  vx = rep(0,n)
  vy = rep(0,n)
  ax = rep(0,n)
  ay = rep(0,n)
  
  x[1] = x0
  y[1] = y0
  z[1] = f(x[1], y[1])
  t[1] = t0
  vx[1] = vx0
  vy[1] = vy0
  ax[1] = 0
  ay[1] = 0
  
  
  for (i in 2:n) {
    t[i] = t[i-1] + dt
    R = ( ( f.dash(x[i-1], y[i-1]) * (vx[i-1]^2 + vy[i-1]^2 - ((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2 ) / u(x[i-1],y[i-1]) ) + (((x[i-1]*vx[i-1]+y[i-1]*vy[i-1])/u(x[i-1],y[i-1]))^2)*(-1/(4*(u(x[i-1],y[i-1])-1)^1.5)) + g)/((u(x[i-1],y[i-1]))*(f.dash(x[i-1], y[i-1]) + (1/f.dash(x[i-1], y[i-1])))) 
    ax[i] = -x[i-1]*R
    ay[i] = -y[i-1]*R
    vx[i] = vx[i-1] + ax[i-1]*dt 
    vy[i] = vy[i-1] + ay[i-1]*dt
    x[i] = x[i-1] + vx[i-1]*dt + ax[i-1]*(dt^2)/2
    y[i] = y[i-1] + vy[i-1]*dt + ay[i-1]*(dt^2)/2
    z[i] = f(x[i], y[i])
    
    # Sys.sleep(0.01)
  }
  
  for (i in seq(1, n, by = 50)) {
    p <- plot_ly(source = "plot2") %>%
      add_surface(x = ~X, y = ~Y, z = ~Z) %>%
      add_markers(x = ~x[1:i], y = ~y[1:i], z = ~z[1:i], opacity = 0.8, color = c("red"), size = 2)
    orca(p, paste0("p",i,".png"))
  }
  
}
gravity.well(0, 10, 0, 0, 5, 5000, 0.01)

library(gifski)
png_files <- list.files(getwd(), pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 0.1)
