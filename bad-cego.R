M = read.csv("../instance_0.9",header = FALSE)
M
n=10
perm= c(1:10)

perm = PerMallows::runif.permutation(1,n) 
inv = inverse.perm(perm)
inv
lop_fitness<-function(perm,mat){
  n = length(perm)
  inve = inverse.perm(perm)
  sol = 0
  for(i in 1:n){
    for(j in 1:n){
      if (i<j){#odio r
        sol = sol+ M[inve[i],inve[j]]
      }
    }
  }
  return(sol);
}
perm
lop_fitness(perm,M)
