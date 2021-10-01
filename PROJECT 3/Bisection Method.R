# BISECTION METHOD

l = 0
r = pi/2
f =function(x) cos(x) - x
for (i in 1:20) {
  m = (l+r)/2
  if (f(m)*f(l) < 0) {
    r = m
  }
  else{
    l = m
  }
  cat(l, r, '\n')
}


# EXERCISE 1
bis = function(f, l, r, n){
  stopifnot(f(l)*f(r) < 0)
  for (i in 1:n) {
    m = (l+r)/2
    if(f(l)*f(m) < 0){
      r = m
    }
    else{
      l = m
    }
    cat(l, r, "\n")
  }
}

# EXERCISE 2
f = function(x) 2*exp(x) - 2*x - 3
bis(f, 0, 2, 50)


# BISECTION + NR COMBO
# EQN: 2*e^x - 2*x - 3 = 0 in (0, 2)
d = function(x) 2*exp(x) - 2
NR = function(f, d, x0, n){
  x = x0
  for (i in 1:n) {
    (x = x - f(x)/d(x))
    print(x)
  }
}

bis.NR = function(f, l, r){
  stopifnot(f(l)*f(r) < 0)
  for (i in 1:5) {
    m = (l+r)/2
    if(f(l)*f(m) < 0){
      r = m
    }
    else{
      l = m
    }
    cat(l, r, "\n")
  }
  NR(f, d, m, 5)
}

