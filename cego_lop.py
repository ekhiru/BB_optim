from rpy2.robjects.packages import STAP
import numpy as np
import itertools as it
import mallows_kendall as mk
import os
os.environ['RPY2_CFFI_MODE'] = "API" # bug in cffi 1.13.0 https://bitbucket.org/rpy2/rpy2/issues/591/runtimeerror-found-a-situation-in-which-we

from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri
numpy2ri.activate()
import rpy2.rinterface as ri

# funcion de distancia entre permutaciones
@ri.rternalize
def kendallTau(A, B):
    n = len(A)
    pairs = it.combinations(range(n), 2)
    distance = 0
    for x, y in pairs:
        a = A[x] - A[y]
        try:
            b = B[x] - B[y]# if discordant (different signs)
        except:
            print("ERROR kendallTau, check b",A, B, x, y)
        if (a * b < 0):
            distance += 1
    return distance

def synthetic_LOP(n, m, phi):
  instance = np.zeros((n,n))
  s = np.array(mk.samplingMM(m,n, phi=phi, k=None))
  for i in range(n):
      for j in range(i+1,n):
          instance[i,j] = (s[:,i]< s[:,j]).sum()
          instance[j,i] = m - instance[i,j]
  return instance

def u_phi(sample,s0, ws):
    m , n = np.array(sample).shape
    #if s0 is None: s0 = np.argsort(np.argsort(rankings.sum(axis=0))) #borda
    dist_avg = np.array([mk.kendallTau(perm, s0) for perm in sample]*ws).sum()/ws.sum() #np.mean(np.array([kendallTau(s0, perm) for perm in rankings]))
    try:
        theta = optimize.newton(mk.mle_theta_mm_f, 0.01, fprime=mk.mle_theta_mm_fdev, args=(n, dist_avg), tol=1.48e-08, maxiter=500, fprime2=None)
    except:
        #if dist_avg == 0.0: return s0, np.exp(-5)#=phi
        print("error. fit_mm. dist_avg=",dist_avg, dist_avg == 0.0)
        print(s0)
        raise
    if theta < 0:
        theta = 0.001
    return np.exp(-theta) # theta = - np.log(phi)

def get_fitness(perm, instance, problem):
    sol = 0
    n = len(perm)#sum inthe upper triangle. we have to maximize this
    inverse = np.argsort(perm)
    for i in range(n):
        for j in range(i,n):
            sol += instance[int(inverse[i]),int(inverse[j])]
    return sol

def uborda(sample,ws):
    mul = (sample*ws[:, None]).sum(axis=0)
    return np.argsort(np.argsort(mul))

def paint(df):
    for problem in df.problem.drop_duplicates():
      sns.set_style("whitegrid")
      color_variable = 'rho'
      y_variables = ['fitnesses','phi_estim','phi_sample','dist']
      for y_variable in y_variables:
        plt.figure(figsize=(15,5))
        palette = sns.color_palette("husl", len(df[color_variable].drop_duplicates()))
        sns.lineplot(x='m',hue=color_variable, y=y_variable, data=df[df.problem==problem], palette=palette)
        plt.show()

def get_expected_distance(iter_ratio, ini_dist):
    return ini_dist-ini_dist*iter_ratio
#get_expected_distance(0, 4)

class LOP:
  def __init__(self, n, m, phi):
    self.instance = synthetic_LOP(n, m, phi)

  # Minimized
  def fitness(self, perm):
    print(perm)
    x = get_fitness(perm, self.instance,"LOP")
    print(f'x = {x}')
    return x

def runR(lop):
    np.random.seed(0)
    #cego = importr("CEGO")
    rstring = """
    library(CEGO)
    options(error=recover)
my_optimCEGO <- function (x = NULL, fun, control = list()) 
{
    con <- list(evalInit = 2, vectorized = FALSE, verbosity = 0, 
        plotting = FALSE, targetY = -Inf, budget = 100, creationRetries = 100, 
        distanceFunction = distancePermutationHamming, creationFunction = solutionFunctionGeneratorPermutation(6), 
        infill = infillExpectedImprovement, model = modelKriging, 
        modelSettings = list(), optimizer = optimEA, optimizerSettings = list(), 
        initialDesign = designMaxMinDist, archiveModelInfo = NULL, 
        initialDesignSettings = list())
    con[names(control)] <- control
    control <- con
    rm(con)
    count <- control$evalInit
    archiveModelInfo <- control$archiveModelInfo
    vectorized <- control$vectorized
    verbosity <- control$verbosity
    plotting <- control$plotting
    creationFunction <- control$creationFunction
    distanceFunction <- control$distanceFunction
    if (is.null(control$initialDesignSettings$distanceFunction)) 
        control$initialDesignSettings$distanceFunction <- distanceFunction
    fun
    if (!vectorized) 
        fn <- function(x) {
    y <- lapply(x, fun)
    print(y)
    return(unlist(y))
    }
    else fn <- fun
    res <- list(xbest = NA, ybest = NA, x = NA, y = NA, distances = NA, 
        modelArchive = NA, count = count, convergence = 0, message = "")
    msg <- "Termination message:"
    res$x <- control$initialDesign(x, creationFunction, count, 
        control$initialDesignSettings)
    distanceHasParam <- FALSE
    if (is.function(distanceFunction)) {
        if (length(distanceFunction) == 1) 
            distanceHasParam <- length(formalArgs(distanceFunction)) > 
                2
        else distanceHasParam <- any(sapply(sapply(distanceFunction, 
            formalArgs, simplify = FALSE), length) > 2)
        if (!distanceHasParam) 
            res$distances <- distanceMatrixWrapper(res$x, distanceFunction)
    }
    print(fun)
    print(fn)
    print(res$x)
    res$y <- fn(res$x)
cat("res$y: ")
print(res$y)
    indbest <- which.min(res$y)
cat("indbest: ", indbest, "\n")
    res$ybest <- res$y[[indbest]]
    res$xbest <- res$x[[indbest]]
    model <- buildModel(res, distanceFunction, control)
    if (!is.null(archiveModelInfo)) {
        res$modelArchive <- list()
        archiveIndex <- 1
        if (identical(model, NA)) {
            res$modelArchive[[archiveIndex]] <- rep(NA, length(archiveModelInfo))
            names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
        }
        else {
            res$modelArchive[[archiveIndex]] <- model$fit[archiveModelInfo]
            names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
        }
    }
    useEI <- is.function(control$infill)
    while ((res$count < control$budget) & (res$ybest > control$targetY)) {
        if (!identical(model, NA)) {
            optimres <- CEGO:::optimizeModel(res, creationFunction, 
                model, control)
            duplicate <- list(optimres$xbest) %in% res$x
            improved <- optimres$ybest < optimres$fpredbestKnownY
        }
        else {
            msg <- paste(msg, "Model building failed, optimization stopped prematurely.")
            warning("Model building failed in optimCEGO, optimization stopped prematurely.")
            res$convergence <- -1
            break
        }
        res$count <- res$count + 1
        if (!duplicate && ((improved || useEI))) {
            res$x[[res$count]] <- optimres$xbest
        }
        else {
            if (!distanceHasParam) {
                designSize <- length(res$x) + 1
                if (is.list(distanceFunction)) 
                  dfun <- distanceFunction[[1]]
                else dfun <- distanceFunction
                xc <- designMaxMinDist(res$x, creationFunction, 
                  designSize, control = list(budget = control$creationRetries, 
                    distanceFunction = dfun))
                res$x[[res$count]] <- xc[[designSize]]
            }
            else {
                res$x[[res$count]] <- optimres$xbest
            }
        }
        res$x <- removeDuplicates(res$x, creationFunction)
        res$y <- c(res$y, fn(res$x[res$count]))
        indbest <- which.min(res$y)
        res$ybest <- res$y[[indbest]]
        res$xbest <- res$x[[indbest]]
        if (verbosity > 0) 
            print(paste("Evaluations:", res$count, "    Quality:", 
                res$ybest))
        if (plotting) {
            plot(res$y, type = "l", xlab = "number of evaluations", 
                ylab = "y")
            abline(res$ybest, 0, lty = 2)
        }
        if (!distanceHasParam & is.function(distanceFunction)) 
            res$distances <- distanceMatrixUpdate(res$distances, 
                res$x, distanceFunction)
        model <- buildModel(res, distanceFunction, control)
        if (!is.null(archiveModelInfo)) {
            archiveIndex <- archiveIndex + 1
            if (identical(model, NA)) {
                res$modelArchive[[archiveIndex]] <- rep(NA, length(archiveModelInfo))
                names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
            }
            else {
                res$modelArchive[[archiveIndex]] <- model$fit[archiveModelInfo]
                names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
            }
        }
    }
    if (min(res$ybest, na.rm = TRUE) <= control$targetY) {
        msg <- paste(msg, "Successfully achieved target fitness.")
        res$convergence <- 1
    }
    else if (res$count >= control$budget) {
        msg <- paste(msg, "Target function evaluation budget depleted.")
    }
    res$message <- msg
    res$distances <- NULL
    res
}
    my_cego <- function(fun, dist, budget = 15)
    {
    seed <- 0
    #distance
    #dF <- distancePermutationHamming
    #mutation
    mF <- mutationPermutationSwap
    #recombination
    rF <- recombinationPermutationCycleCrossover
    #creation
    cF <- function()sample(6)
    #start optimization
    set.seed(seed)
    res1 <- my_optimCEGO(x = NULL,
                      fun = fun,
                      control = list(creationFunction=cF,
                                     distanceFunction = dist,
                                     optimizerSettings=list(budget=100,popsize=10,
                                                            mutationFunction=mF,
                                                            recombinationFunction=rF),
                                     evalInit=5,budget=budget,targetY=0,verbosity=1,
                                     model=modelKriging,
                                     vectorized=FALSE))
    print(res1)
    return(list(res1$xbest, res1$ybest))
    }
    """

    # with open('myfunc.r', 'r') as f:
    #     rstring = f.read()
    rcode = STAP(rstring, "rcode")
    best_x, best_fitness = rcode.my_cego(ri.rternalize(r_lop_fitness),
                                         dist = kendallTau)
    best_x = np.asarray(best_x)
    best_fitness = np.asarray(best_fitness)[0]
    print(f'best: {best_x}\nbest_fitness: {best_fitness}')

    
lop = LOP(6,100, phi=0.9)
@ri.rternalize
def r_lop_fitness(x):
    y = lop.fitness(x)
    return ri.FloatSexpVector([y])


y = runR(lop)
