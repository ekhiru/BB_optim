library(PerMallows)
library(CEGO)

distKendall <- function(x, y)
{
  n <- length(x)
  stopifnot(length(x)== length(y))
  pairs <- t(combn(1:n,2))
  a <- x[pairs[,1]] - x[pairs[,2]] 
  b <- y[pairs[,1]] - y[pairs[,2]]
  sum( a * b < 0)
}

gen_lop_fitness <- function(mat) {
  force(mat)
  n <- nrow(mat)
  lop_fitness <- function(p) {
    z <- 0
    for (i in 1:(n-1)) {
      z <- z + sum(mat[p[i], p[(i+1):n]])
    }
    stopifnot(!is.na(z))
    #print(z)
    -z
  }
  lop_fitness
}

my_cego <- function(fun, dist, n, m_ini = 5, budget = 400, seed = 0, budgetGA = 10000)
{
  set.seed(seed)
  # mutation
  mF <- mutationPermutationInterchange # changed as paper mutationPermutationSwap
  # recombination
  rF <- recombinationPermutationCycleCrossover
  #creation
  cF <- function() sample(n)
  # start optimization
  res <- optimCEGO(x = NULL,
                   fun = fun,
                   control = list(creationFunction=cF,
                                  #distanceFunction = dist,
                                  optimizerSettings=list(budget=budgetGA,popsize=20,
                                                         mutationFunction=mF,
                                                         recombinationFunction=rF),
                                  evalInit=m_ini,budget=budget,verbosity=1,
                                  model=modelKriging,
                                  vectorized=FALSE))
  print(res)
  return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
}


instance <- "instance_0.9"
M <- read.csv(instance,header = FALSE)
n <- nrow(M)
perm <- sample(1:n)

lop_fitness <- gen_lop_fitness(M)
print(lop_fitness(1:n))
print(lop_fitness(n:1))

my_cego(fun = lop_fitness, dist = distKendall, #PerMallows::distance, # Kendall distance
        n = n)

