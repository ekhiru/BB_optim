import os
#os.environ['RPY2_CFFI_MODE'] = "API" # bug in cffi 1.13.0 https://bitbucket.org/rpy2/rpy2/issues/591/runtimeerror-found-a-situation-in-which-we

import numpy as np
from mallows_kendall import kendallTau
import pandas as pd

from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import STAP
numpy2ri.activate()
import rpy2.rinterface as ri

# funcion de distancia entre permutaciones
@ri.rternalize
def r_kendallTau(A,B):
    return kendallTau(A,B)

def cego(instance, seed, budget, # m: number of evaluations
         m_ini, budgetGA):
    # Reset the list of recorded evaluations.
    instance.reset()
    
    rstring = """
    library(CEGO)
    # This is identical to the function in the package but takes also a maxTime parameter.
optimCEGO <- function (x = NULL, fun, control = list()) 
{
    con <- list(evalInit = 2, vectorized = FALSE, verbosity = 0, 
        plotting = FALSE, targetY = -Inf, budget = 100, creationRetries = 100, 
        distanceFunction = distancePermutationHamming, creationFunction = solutionFunctionGeneratorPermutation(6), 
        infill = infillExpectedImprovement, model = modelKriging, 
        modelSettings = list(), optimizer = optimEA, optimizerSettings = list(), 
        initialDesign = designMaxMinDist, archiveModelInfo = NULL, 
        initialDesignSettings = list(), maxTime = 3600 * 24 * 3)
    con[names(control)] <- control
    control <- con
    rm(con)
    maxTime <- control$maxTime + proc.time()["elapsed"] 
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
        fn <- function(x) unlist(lapply(x, fun))
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
            res$distances <- CEGO:::distanceMatrixWrapper(res$x, distanceFunction)
    }
    res$y <- fn(res$x)
    indbest <- which.min(res$y)
    res$ybest <- res$y[[indbest]]
    res$xbest <- res$x[[indbest]]
    model <- CEGO:::buildModel(res, distanceFunction, control)
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
    while ((res$count < control$budget) & (res$ybest > control$targetY) & ((maxTime - proc.time()["elapsed"]) > 0)) {
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
                xc <- CEGO:::designMaxMinDist(res$x, creationFunction, 
                  designSize, control = list(budget = control$creationRetries, 
                    distanceFunction = dfun))
                res$x[[res$count]] <- xc[[designSize]]
            }
            else {
                res$x[[res$count]] <- optimres$xbest
            }
        }
        res$x <- CEGO::removeDuplicates(res$x, creationFunction)
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
            res$distances <- CEGO:::distanceMatrixUpdate(res$distances, 
                res$x, distanceFunction)
        model <- CEGO:::buildModel(res, distanceFunction, control)
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
    else if ((maxTime - proc.time()["elapsed"]) <= 0) {
        msg <- paste(msg, "maxTime reached.")
    }
    res$message <- msg
    res$distances <- NULL
    res
    }

    my_cego <- function(fun, dist, n, m_ini = 5, budget = 15, seed = 0, budgetGA = 100)
    {
    set.seed(seed)
    # mutation
    mF <- mutationPermutationInterchange # changed as paper mutationPermutationSwap
    # recombination
    rF <- recombinationPermutationCycleCrossover
    #creation
    cF <- function() sample(n)
    # start optimization
#    print("antes del optimCEGO")
    res <- optimCEGO(x = NULL,
                     fun = fun,
                     control = list(creationFunction=cF,
                                     distanceFunction = dist,
                                     #distanceFunction = distancePermutationSwap,
                                     optimizerSettings=list(budget=budgetGA,popsize=20,
                                                            mutationFunction=mF,
                                                            recombinationFunction=rF),
                                     evalInit=m_ini,budget=budget,verbosity=1,
                                     model=modelKriging,
                                     vectorized=FALSE))

    return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
    }
    """
    rcode = STAP(rstring, "rcode")
    r_fitness = instance.make_r_fitness()
    best_x, best_fitness, x, y = rcode.my_cego(r_fitness,
                                                dist = r_kendallTau,
                                                n = instance.n,
                                                m_ini = m_ini,
                                                budget = budget,
                                                seed = seed,
                                                budgetGA = budgetGA)

    # We use what instance recorded because CEGO may not get us what was
    # actually evaluated.
    return pd.DataFrame(dict(
        Fitness=instance.evaluations, x=instance.solutions, m_ini = m_ini, seed = seed, budget = budget, budgetGA = budgetGA,
        Distance = [ instance.distance_to_best(perm) for perm in instance.solutions]))
    

# rep : replication number.
def runCEGO(instance, m_ini, m, rep, best_known_sol, worst_known_sol, budgetGA):
    
    rstring = """
    library(CEGO)

    my_cego <- function(fun, dist, n, m_ini = 5, budget = 15, seed = 0, budgetGA = 100)
    {
    set.seed(seed)
    # mutation
    mF <- mutationPermutationInterchange # changed as paper mutationPermutationSwap
    # recombination
    rF <- recombinationPermutationCycleCrossover
    #creation
    cF <- function() sample(n)
    # start optimization
#    print("antes del optimCEGO")
    res <- optimCEGO(x = NULL,
                     fun = fun,
                     control = list(creationFunction=cF,
                                     distanceFunction = dist,
                                     #distanceFunction = distancePermutationSwap,
                                     optimizerSettings=list(budget=budgetGA,popsize=20,
                                                            mutationFunction=mF,
                                                            recombinationFunction=rF),
                                     evalInit=m_ini,budget=budget,verbosity=1,
                                     model=modelKriging,
                                     vectorized=FALSE))

    return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
    }
    """

    # with open('myfunc.r', 'r') as f:
    #     rstring = f.read()
    rcode = STAP(rstring, "rcode")
    r_fitness = instance.make_r_fitness()
    best_x, best_fitness, x , y = rcode.my_cego(r_fitness,
                                                dist = r_kendallTau,
                                                n = instance.n,
                                                m_ini = m_ini,
                                                budget = m,
                                                seed = rep,
                                                budgetGA = budgetGA)

    x = np.asarray(x)
    y = np.asarray(y)
    #best_x = np.asarray(best_x)
    #best_fitness = np.asarray(best_fitness)[0]
    #print(f'best: {best_x}\nbest_fitness: {best_fitness}')
    df = pd.DataFrame()#columns=['problem','rep','m','rho','fitnesses','phi_estim','phi_sample','dist'])
    # Evaluate without saving the solution.
    if worst_known_sol is None:
        worst_known_fit = 0
    else:
        worst_known_fit = instance.fitness_nosave(worst_known_sol)
    
    # FIXME: We should move all these calculations outside here so that we don't need the best known or worst-known here.
    if best_known_sol is None:
        df['Fitness'] = instance.evaluations
    else:
        best_known_fit = instance.fitness_nosave(best_known_sol)
        df['Fitness'] = (instance.evaluations - best_known_fit) / np.abs(worst_known_fit - best_known_fit)
    df['Problem'] = instance.problem_name
    df['Solver'] = 'CEGO'
    df['Sample size'] = range(m)
    df['rep'] = rep
    #    df['budgetGA'] = budgetGA this must be set for all, including uMM, so that we can filter appropriately
    # FIXME: We should calculate the distance outside here so we do not need the best_known_sol here.
    if best_known_sol is not None:
        df['Distance'] = [ kendallTau(perm, best_known_sol) / (instance.n * (instance.n - 1) / 2) for perm in instance.solutions]
    return df
