classicga <- function(fitness,
                      min,
                      max,
                      populationSize = 50,
                      crossoverProb = 0.2,
                      mutationProb = 0.1,
                      iterations = 100,
                      elitismPercentage = 0.1,
                      verbose = FALSE,
                      convergenceIters = Inf,
                      initialPopulation = NULL
) {
  # validate arguments
  if(missing(fitness) | !is.function(fitness)) {
    stop("A fitness function is required")
  }
  if(missing(min) | missing(max)) {
    stop("The min and max vectors are required")
  }
  if(iterations < 1) {
    stop("The number of iterations must be greater than 0")
  }
  if(crossoverProb < 0 | crossoverProb > 1) {
    stop("The provided crossover probability is not valid")
  }
  if(mutationProb < 0 | mutationProb > 1) {
    stop("The provided mutation probability is not valid")
  }
  if(elitismPercentage < 0 | elitismPercentage > 1) {
    stop("The provided elitism percentage is not valid")
  }
  if(length(min) != length(max)) {
    stop("The min and max vectors must have the same length")
  }
  
  dimension <- length(min)
  fitnessVec <- rep(NA, populationSize)
  eps = sqrt(.Machine$double.eps)
  elitism <- base::max(1, round(populationSize * elitismPercentage))
  fitnessSummary <- matrix(as.double(NA), nrow = iterations, ncol = 6)
  run <- 1
  
  
  if(is.null(initialPopulation))
  { initialPopulation <- matrix(nrow = 0, ncol = dimension) }
  else
  { if(is.vector(initialPopulation)) 
  { if(dimension > 1) initialPopulation <- matrix(initialPopulation, nrow = 1)
  else          initialPopulation <- matrix(initialPopulation, ncol = 1) }
    else
    { initialPopulation <- as.matrix(initialPopulation) }
    if(dimension != ncol(initialPopulation))
      stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
  }
  
  population <- matrix(as.double(NA), nrow = populationSize, ncol = dimension)

  ng <- min(nrow(initialPopulation), populationSize)
  if(ng > 0) { 
    population[1:ng,] <- initialPopulation 
  }
  if(populationSize > ng) { 
    population[(ng+1):populationSize,] <- generatePopulation(min, max, populationSize)[1:(populationSize-ng),]
  }

  if(verbose) print(population)
  
  for (iter in seq_len(iterations)) {
    # evaluate fitness
    for(i in seq_len(populationSize)) {
      if(is.na(fitnessVec[i])) {
        fitnessVec[i] <- fitness(population[i,])
      }
    }
    
    # update summary
    x <- na.exclude(as.vector(fitnessVec))
    summaryFiveNum <- fivenum(x)
    fitnessSummary[iter,] <- c(max = summaryFiveNum[5], mean = mean(x), q3 = summaryFiveNum[4], median = summaryFiveNum[3], q1 = summaryFiveNum[2], min = summaryFiveNum[1])
    if(iter > 1) {
      run <- compareConvergence(fitnessSummary[seq(iter),1])
    }
    if(run >= convergenceIters) break
    #if(max(fitnessVec, na.rm = TRUE) >= Inf) break
    
    
    if(iter == iterations) break
    
    ord <- order(fitnessVec, decreasing = TRUE)
    populationSorted <- population[ord,,drop=FALSE]
    fitnessSorted <- fitnessVec[ord]

    # selection
    selected <- rouletteSelection(fitnessVec, population, populationSize)
    population <- selected$population
    fitnessVec <- selected$fitness
    if(verbose) print(selected)
    
    # crossover
    if(crossoverProb > 0) {
      nmating <- floor(populationSize / 2)
      mating <- matrix(sample(1:populationSize, size = (2 * nmating)), ncol = 2)
      for(i in seq_len(nmating)) { 
        if(crossoverProb > runif(1)) { 
          parents <- mating[i,]
          crossoverResult <- singlePointCrossover(fitnessVec, population, parents)
          population[parents,] <- crossoverResult$children
          fitnessVec[parents] <- crossoverResult$fitness
        }
      }
    }
    if(verbose) print(population)

    # mutation
    if(mutationProb > 0) {
      for(i in seq_len(populationSize)) {
        mutated <- gaussianMutation(population[i,], mutationProb, min, max)
        # if(!identical(mutated, population[i])) {
          fitnessVec[i] <- NA
        # }
        population[i,] <- mutated
      }
    }
    if(verbose) print(population)
    
    # elitism
    ord <- order(fitnessVec, na.last = TRUE)
    u <- which(!duplicated(populationSorted, margin = 1))
    population[ord[1:elitism],] <- populationSorted[u[1:elitism],]
    fitnessVec[ord[1:elitism]] <- fitnessSorted[u[1:elitism]]
  }
  
  fitnessValue <- max(fitnessVec, na.rm = TRUE)
  bestResult <- which(fitnessVec == fitnessValue)
  solution <- population[bestResult,,drop=FALSE]
  return(list(solution=solution,population=population,fitness=fitnessVec))
}

rouletteSelection <- function(fitness, population, popSize) {
  prob <- abs(fitness)/sum(abs(fitness))
  sel <- sample(1:popSize, size = popSize, 
                prob = pmin(pmax(0, prob), 1, na.rm = TRUE),
                replace = TRUE)
  result <- list(population = population[sel,,drop=FALSE],
              fitness = fitness[sel])
  return(result)
}

singlePointCrossover <- function(fitnessVec, population, parents) {
  fitness <- fitnessVec[parents]
  parents <- population[parents,,drop = FALSE]
  n <- ncol(parents)
  
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  fitnessChildren <- rep(NA, 2)
  crossOverPoint <- sample(0:n, size = 1)
  
  if (crossOverPoint == 0) { 
    children[,] <- parents[2:1,]
    fitnessChildren <- rev(fitness)
  } else if (crossOverPoint == n) { 
    children <- parents
    fitnessChildren <- fitness 
  } else { 
    children[1,] <- c(parents[1,1:crossOverPoint],
                      parents[2,(crossOverPoint + 1):n])
    children[2,] <- c(parents[2,1:crossOverPoint],
                    parents[1,(crossOverPoint + 1):n])
  }
  result <- list(children = children, fitness = fitnessChildren)
  return(result)
}

gaussianMutation <- function(solution, prob, min, max) {
  n <- length(solution)
  mutant <- solution
  # select variables for mutation
  idx <- which(runif(n) < prob)
  mutate <- rnorm(length(idx), mean = 0, sd = 0.04)
  mutant[idx] <- solution[idx] + mutate
  # correct bounds
  mutant <- pmax(pmin(mutant, max), min)
  return(mutant)
}

compareConvergence <- function(x) {
  x <- as.vector(x)
  sum(rev(x) >= (max(x, na.rm = TRUE) - sqrt(.Machine$double.eps)))
}

generatePopulation <- function(min, max, populationSize){
  nvars <- length(min)
  population <- matrix(as.double(NA), nrow = populationSize, ncol = nvars)
  for(j in 1:nvars){ 
    population[,j] <- runif(populationSize, min[j], max[j]) 
  }
  return(population)
}