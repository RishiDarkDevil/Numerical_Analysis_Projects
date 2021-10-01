# MAXIMUM LIKELIHOOD ESTIMATOR RABBIT

L = function(p) (p^46)*((1-p)^2 + 2*p*(1-p))^77

L.dash = function(p) (46*(p^45)*((1-p)^2 + 2*p*(1-p))^77) + ((p^46)*(77*((1-p)^2 + 2*p*(1-p))^76)*(-2*p))

L.ddash = function(p) (46*45*L(p)/p) + ((-2*46*77)*L(p)/(1-p^2)) + ((-2*77*47)*L(p)/(1-p^2)) + ((2*2*77*76)*(p^2)*(L(p)/((1-p^2)^2)))

NR = function(f, d, x0, n){
  x = x0
  for (i in 1:n) {
    (x = x - f(x)/d(x))
    print(x)
  }
}

NR(L.dash, L.ddash, 0.5, 1000)

goldsectmax <- function (f, a, b, tol = 1e-3, m = 100) {
  iter <- 0
  phi <- ( sqrt (5) - 1) / 2
  a.star <- b - phi * abs(b - a)
  b.star <- a + phi * abs(b - a)
  while (abs(b - a) > tol) {
    iter <- iter + 1
    if ( iter > m) {
      warning (" iterations maximum exceeded ")
      break
    }
    if(f(a.star ) > f(b.star )) {
      b <- b.star
      b.star <- a.star
      a.star <- b - phi * abs(b - a)
    } else {
      a <- a.star
      a.star <- b.star
      b.star <- a + phi * abs(b - a)
    }
  }
  return ((a + b) / 2)
}

goldsectmax(L, 0, 1, 0.0000001, 34)


# Visualizing NR Iterations
NR = function(f, d, x0, n, tol){
  p1 <- function_plots %>%
    ggplot(aes(p, L.values)) +
    geom_line(size = 2, color = "blue") +
    geom_hline(yintercept = 0) +
    ggtitle("Likelihood Function") +
    theme_bw() +
    scale_y_continuous(labels = NULL)
  p2 <- function_plots %>%
    ggplot(aes(p, L.dash.values)) +
    geom_line(size = 2, color = "blue") +
    geom_hline(yintercept = 0) +
    ggtitle("Derivative of Likelihood Function") +
    theme_bw() +
    scale_y_continuous(labels = NULL)
  x = x0
  for (i in 1:n) {
    x_old <- x
    (x = x - f(x)/d(x))
    p2 <- p2 + geom_point(data = tibble(p = x_old, L.dash.val = L.dash(x_old)), aes(p, L.dash.val), size = 1, color = "red") + geom_line(data = tibble(p = c(x_old, x), L.dash.val = c(L.dash(x_old), L.dash(x))), aes(p, L.dash.val), size = 1, color = "orange")
    p1 <- p1 + geom_point(data = tibble(p = x_old, L.val = L(x_old)), aes(p, L.val), size = 1, color = "red") + geom_line(data = tibble(p = c(x_old, x), L.val = c(L(x_old), L(x))), aes(p, L.val), size = 1, color = "orange")
    plt <- annotate_figure(ggarrange(p1, p2, ncol = 2), top = text_grob("NEWTON-RAPHSON'S ITERATIONS", face = "bold", size = 16))
    ggsave(paste0("p",i,".png"), plot = plt, width = 12, height = 0.618*12)
    if(abs(x-x_old) < tol){
      print(p1)
      return(x) 
    }
  }
  print(p1)
  warning("Maximum Iteration Limit Exceeded!")
  return(x)
}

NR(L.dash, L.ddash, 0.5, 100, 1e-5)

library(gifski)
png_files <- list.files(getwd(), pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 0.2)
