# NEWTON RAPHSON'S METHOD

# EQN: e^x = cosx, for x in (-pi/2, 0)
x = -pi/4
(x = x - (exp(x) - cos(x))/(exp(x) + sin(x)))

# EQN: 1 - 2x + 3x^3 + 5x^4 - x^5 = 0
x = 1
(x = x - (1 - 2*x + 3*(x^3) + 5*(x^4) - (x^5)) / (-2 + (9*x^2) + (20*x^3) - (5*x^4)))

NR = function(f, d, x0, n){
  x = x0
  (x = x - f(x)/d(x))
}

# EXERCISE 2

# Example function- EQN: 1 - 2x + 3x^3 + 5x^4 - x^5 = 0
f = function(x) {1 - 2*x + 3*(x^3) + 5*(x^4) - (x^5)}
f.dash = function(x) {-2 + (9*x^2) + (20*x^3) - (5*x^4)}
start.x = 5
iter = 100

# Function- NR formula: x_i = x_(i-1) - f(x_(i-1)) / f'(x_(i-1))
# f here in the example taken is 1 - 2x + 3x^3 + 5x^4 - x^5 and f' is -2 + (9*x^2) + (20*x^3) - (5*x^4)
NR2 = function(f, d, x0, n){
  for (i in 1:n) {
    x0 = x0 - f(x0)/d(x0)
    print(x0)
  } 
}

# Running the NR2 function with example data
NR2(f, f.dash, start.x, iter)
