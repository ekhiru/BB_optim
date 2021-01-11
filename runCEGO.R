argv <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(argv[1])
eval_ranks <- as.numeric(argv[2])
elapsed <- proc.time()
#seed <- 12345
#eval_ranks <- 1
library(CEGO)
cat("seed = ", seed, "  eval_ranks = ", eval_ranks, "\n")
print(sessionInfo())
# This is identical to the function in the package but takes also a maxTime parameter.
my_optimCEGO <- function (x = NULL, fun, control = list()) 
{
    con <- list(evalInit = 2, vectorized = FALSE, verbosity = 0, 
        plotting = FALSE, targetY = -Inf, budget = 100, creationRetries = 100, 
        distanceFunction = distancePermutationHamming, creationFunction = solutionFunctionGeneratorPermutation(6), 
        infill = infillExpectedImprovement, model = modelKriging, 
        modelSettings = list(), optimizer = optimEA, optimizerSettings = list(), 
        initialDesign = designMaxMinDist, archiveModelInfo = NULL, 
        initialDesignSettings = list(), maxTime = 3600 * 24 * 6, eval_ranks = FALSE)
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
    if (!vectorized) { 
        fn <- if (control$eval_ranks) function(x) unlist(lapply(x, fun)) 
              else function(x) unlist(lapply(x, function(y) fun(order(y) - 1))) 
    } else {
     # fn <- fun
     stop("We do not handle vectorized functions")
    }
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
        if ((!duplicate && (improved || useEI)) || distanceHasParam) {
            res$x[[res$count]] <- optimres$xbest
        } else { # !distanceHasParam
                designSize <- length(res$x) + 1
                if (is.list(distanceFunction)) 
                  dfun <- distanceFunction[[1]]
                else dfun <- distanceFunction
                xc <- CEGO:::designMaxMinDist(res$x, creationFunction, 
                  designSize, control = list(budget = control$creationRetries, 
                    distanceFunction = dfun))
                res$x[[res$count]] <- xc[[designSize]]
                # Modified to print
                print(paste0("Random solution:", res$count))
        }
        res$x <- CEGO::removeDuplicates(res$x, creationFunction)
        res$y <- c(res$y, fn(res$x[res$count]))
        indbest <- which.min(res$y)
        res$ybest <- res$y[[indbest]]
        res$xbest <- res$x[[indbest]]
        if (verbosity > 0) {
            # Modified to print current quality
            print(paste("Evaluations:", res$count, "    Quality:", 
                res$y[res$count], "    Best:", res$ybest))
        }
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

n <- 20
m <- 15
instance <- "rec13"
A <- read.table(text='
78    37   79  98  100  21  82  27    7  22   9  57  84  91   5
13    81  100  77   45  39  60  87   38  91  17  85  43  81  33
47     9   31  40   86  27  69  50   87  34  13  15  95  96  72
26    61    5  22   26  52  57  97   10  68   8  49  41  16  35
68    41   46  58   37  59  22  43   49  21  42  70  13   2  76
47    60   23  23   57  60  35  56   54  73  81  61  15  70  51
63    21   10  50   12  84  91   7   44  75  60  72  28  83  52
52   100   38  79   90  52  81  54   51  11  76  50  24  12  23
91    21   59  55   39  75   3  11   14  24  36  72  48  69  55
00    54   77  53   90  91  53  68   18  86  28  61  86  36  15
63    58   93  24   18  78  74  57  100  92  51  73   6  12  83
96    60   79  84   86  31  41  79   63  17  31  65  41  77  74
54    55   89  77   40  71  21  43   60  58  13   3  73  44  85
50    33   51 100   46  27  30  80   56  50  97  70   9  22  92
17     7   12  66   98  25  76  25   76  25  69  69  73  45  74
46    76   32  32   28  48  63  92   85  69  24  38  57  21  66
91    42   30  39   97  72  96  78   71  75  17  26  94   3  39
99    93   57  23   52  89   4  74   21  10  53  94  59 100  68
20    84   25   1   72  63  79  63   17  38  21  36   1  73  58
18     4   21  43   55  86  99  38   81  74  34  40  61  76 100
')
# create FSP objective function
fun <- benchmarkGeneratorFSP(A,n,m)

set.seed(seed)
# mutation
mF <- mutationPermutationInterchange # changed as paper mutationPermutationSwap
# recombination
rF <- recombinationPermutationCycleCrossover
#creation
cF <- function() sample(n)
# start optimization
#    print("antes del optimCEGO")
budgetGA <- 10^3
budget <- 400
res <- optimCEGO(x = NULL,
                 fun = fun,
                 control = list(creationFunction=cF,
                                distanceFunction = distancePermutationSwap,
                                optimizerSettings=list(budget=budgetGA,popsize=20,
                                                       mutationFunction=mF,
                                                       recombinationFunction=rF),
                                evalInit=10,budget=budget,verbosity=1,
                                model=modelKriging,
                                vectorized=FALSE, eval_ranks = as.logical(eval_ranks)))
elapsed <- proc.time() - elapsed
print(res$y)
print(1:budget)
print(elapsed[3])
# We cannot use budget because optimCEGO may terminate before consuming the budget.
cegores <- data.frame(Fitness=res$y, Instance=instance, Evals = seq(1,length(res$y)), eval_ranks = eval_ranks, seed=seed, time=elapsed[[3]])

saveRDS(cegores, paste0("cegores-er", eval_ranks, "-r", seed, ".rds"))
