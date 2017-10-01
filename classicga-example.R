library("GA")
source("classicga.R")
source("islandsga.R")

# Example 1
f <- function(x)  (x^2+x)*cos(x) # -10 < x < 10
curve(f, -10, 10)
GA <- islandsga(fitness =  f, min = -10, max = 10, populationSize = 500, verbose = FALSE, iterations = 50, 
                islandsCount = 5, hierarchical = FALSE)
print(GA)

GA <- ga(type = "real-valued", fitness = f, min = -10, max = 10, monitor = NULL, maxiter = 1000, run= 50)
summary(GA)
plot(GA)

# Example 2 (difficult function)
Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20)
#filled.contour(x1, x2, f, color.palette = jet.colors)
GA <- islandsga(fitness =  function(x) -Rastrigin(x[1], x[2]),
                min = c(-5.12, -5.12), max = c(5.12, 5.12), populationSize = 500, verbose = FALSE, iterations = 500, 
                islandsCount = 5, hierarchical = FALSE)
#GA <- classicga(fitness =  function(x) -Rastrigin(x[1], x[2]),
#                min = c(-5.12, -5.12), max = c(5.12, 5.12), verbose = FALSE, iterations = 20000)
print(GA)

GA <- GA::ga(type = "real-valued", fitness =  function(x) -Rastrigin(x[1], x[2]),
             min = c(-5.12, -5.12), max = c(5.12, 5.12), 
             popSize = 50, maxiter = 1000)
summary(GA)