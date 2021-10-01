library(rio)
library(rgl)

data = import("data.csv")

n = 996

product = 1
summation = 0
datalist = numeric(996)
dataoutput = numeric(996)
k = 1

L = function(x, p, a) ((a^p)*(x^(p-1))*exp((-1)*a*x))/gamma(p)

for (i in 1:166) {
  for (j in 1:6) {
    product = product*data[i, j]
    summation = summation + data[i, j]
    datalist[k] = data[i, j]
    dataoutput[k] = L(datalist[k],1.946419, 2.878890)
    k = k+1
  }
}

plot(x=datalist, y=dataoutput)

Log.Likelihood <- function(data, p, a) n*log(a^p / gamma(p)) + (p-1)*sum(log(data)) - a*sum(data)

L.dash = function(x){
  p = x[1]
  a = x[2]

  dellogL.delp = n*log(a) + log(product) - n*digamma(p)
  dellogL.dela = (n*p/a) - summation
  
  return(c(dellogL.delp, dellogL.dela))
}

J = function(x) {
  p = x[1]
  a = x[2]

  vec = c(-n*trigamma(p), n/a, n/a, -n*p/a^2)
  J <- matrix(data = vec, ncol = 2)
  return(J)
}

NR2 = function(f, d, x0, n){
  for (i in 1:n) {
    ( x0 = x0 - solve(d(x0), f(x0)) )
    print(x0)
  }  
}

NR2(L.dash, J, c(1,1), 1000)

tibble(x = datalist) %>%
  ggplot(aes(x, ..density..)) +
  geom_histogram(bins = 30) +
  geom_line(data = tibble(x = seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), fit = dgamma(seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), 1.946419, 2.878890)), aes(x, fit), size = 2, color = "blue") +
  theme_bw() +
  ggtitle("Gamma Distribution Fitted")
  

# Visualizing NR Iterations
NR2 = function(f, d, x0, n, tol){
  p1 <- tibble(x = datalist) %>%
    ggplot() +
    geom_histogram(aes(x, ..density..), bins = 30) +
    theme_bw() +
    ggtitle("NR Iterations")
  x = x0
  a <- x[1]
  b <- x[2]
  ggsave(paste0("p","0",".png"), plot = p1 + geom_line(data = tibble(x = seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), fit = dgamma(seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), a, b)), aes(x, fit), size = 2, color = "blue"), width = 12, height = 0.618*12)
  for (i in 1:n) {
    x_old <- x
    ( x = x - solve(d(x), f(x)) )
    print(x)
    a <- x[1]
    b <- x[2]
    ggsave(paste0("p",i,".png"), plot = p1 + geom_line(data = tibble(x = seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), fit = dgamma(seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), a, b)), aes(x, fit), size = 2, color = "blue"), width = 12, height = 0.618*12)
    if(sqrt(sum(x-x_old)^2) < tol){
      print(p1 + geom_line(data = tibble(x = seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), fit = dgamma(seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), a, b)), aes(x, fit), size = 2, color = "blue"))
      return(x) 
    }
  }
  print(p1 + geom_line(data = tibble(x = seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), fit = dgamma(seq(min(datalist, na.rm = TRUE), max(datalist, na.rm = TRUE), 0.01), a, b)), aes(x, fit), size = 2, color = "blue"))
  warning("Maximum Iteration Limit Exceeded!")
  return(x)
}

NR2(L.dash, J, c(1,1), 100, 1e-5)

library(gifski)
png_files <- list.files(getwd(), pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 0.5)
