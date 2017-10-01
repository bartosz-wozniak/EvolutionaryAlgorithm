library("bbob")
source("classicga.R")
source("islandsga.R")

getClassicGaOptimizer <- function(popSize, mutationProb, crossoverProb) {
  myOptimizer <- function(par, fun, lower, upper, max_eval) {
    classicga(fitness=fun, min=lower, max=upper,
              populationSize = popSize, 
              mutationProb = mutationProb,
              crossoverProb = crossoverProb, 
              iterations = max_eval,
              convergenceIters = max_eval,
              elitismPercentage = 0.05)
  }
  return(myOptimizer)
}

getIslandsOptimizer <- function(popSize, mutationProb, crossoverProb, migrationRate, 
                                migrationInterval, islandsCount, hierarchical) {
  myOptimizer <- function(par, fun, lower, upper, max_eval) {
    islandsga(fitness=fun, min=lower, max=upper,
              populationSize = popSize,
              mutationProb = mutationProb,
              crossoverProb = crossoverProb,
              islandsCount = islandsCount,
              migrationInterval = migrationInterval,
              migrationRate = migrationRate,
              iterations = max_eval,
              hierarchical = hierarchical,
              elitismPercentage = 0.05)
  }
}

populationSizes <- c(50, 100, 200)
iterations <- c(100, 500, 1000)
mutationProbs <- c(0.01)
crossoverProbs <- c(0, 0.8)

i <- 1
basicCombinations <- length(populationSizes) * length(iterations) * length(mutationProbs) * length(crossoverProbs)
print("Testing classic GA")
print(paste("Combinations to check:", basicCombinations))
for(popSize in populationSizes) {
  for(itersCount in iterations) {
    for(mutProb in mutationProbs) {
      for(crossoverProb in crossoverProbs) {
        print(paste("Run", i))
        i <- i + 1

        name <- paste0("classicga_", popSize, "_", itersCount, "_", mutProb, "_", crossoverProb)
        print(paste0("Computing ", name))
        print(paste("populationSize", popSize, "iterations", itersCount,
                    "mutationProb", mutProb, "crossOverProb", crossoverProb))
        optimizer <- getClassicGaOptimizer(popSize = popSize, mutationProb = mutProb, crossoverProb = crossoverProb)
        
        oldTime <- Sys.time() 
        bbo_benchmark(optimizer, "classicga", name, budget = itersCount)
        newTime <- Sys.time() - oldTime
        print(newTime)
        cat("\n")
      }
    }
  }
}

migrationRates <- c(0.1, 0.2)
migrationIntervals <- c(10, 50, 100)
islCounts <- c(5, 10)

i <- 1
print("Testing stepping stones island model")
print(paste("Combinations to check:", basicCombinations * length(migrationRates) * length(migrationIntervals) * length(islCounts)))
for(popSize in populationSizes) {
  for(itersCount in iterations) {
    for(mutProb in mutationProbs) {
      for(crossoverProb in crossoverProbs) {
        for(migrationRate in migrationRates) {
          for(migrationInterval in migrationIntervals) {
            for(islandCount in islCounts) {
              print(paste("Run", i))
              i <- i + 1

              name <- paste0("stones_", popSize, "_", itersCount, "_", mutProb, "_", crossoverProb, "_",
                             migrationRate, "_", migrationInterval, "_", islandCount)
              print(paste0("Computing ", name))
              print(paste("populationSize", popSize, "iterations", itersCount,
                          "mutationProb", mutProb, "crossOverProb", crossoverProb, "migrationRate",
                          migrationRate, "migrationInterval", migrationInterval, "islandsCount", islandCount))
              optimizer <- getIslandsOptimizer(popSize = popSize,
                                               mutationProb = mutProb,
                                               crossoverProb = crossoverProb,
                                               migrationRate = migrationRate,
                                               migrationInterval = migrationInterval,
                                               islandsCount = islandCount,
                                               hierarchical = FALSE)
              
              oldTime <- Sys.time() 
              bbo_benchmark(optimizer, "stepping-stones", name, budget = itersCount)
              newTime <- Sys.time() - oldTime
              print(newTime)
              cat("\n")
            }
          }
        }
      }
    }
  }
}

i <- 1
print("Testing hierarchical island model")
print(paste("Combinations to check:", basicCombinations * length(migrationRates) * length(migrationIntervals)))
for(popSize in populationSizes) {
  for(itersCount in iterations) {
    for(mutProb in mutationProbs) {
      for(crossoverProb in crossoverProbs) {
        for(migrationRate in migrationRates) {
          for(migrationInterval in migrationIntervals) {
            print(paste("Run", i))
            i <- i + 1

            name <- paste0("hierarchy_", popSize, "_", itersCount, "_", mutProb, "_", crossoverProb, "_",
                           migrationRate, "_", migrationInterval)
            print(paste0("Computing ", name))
            print(paste("populationSize:", popSize, "iterations", itersCount,
                        "mutationProb", mutProb, "crossOverProb", crossoverProb, "migrationRate",
                        migrationRate, "migrationInterval", migrationInterval))
            optimizer <- getIslandsOptimizer(popSize = popSize,
                                             mutationProb = mutProb,
                                             crossoverProb = crossoverProb,
                                             migrationRate = migrationRate,
                                             migrationInterval = migrationInterval,
                                             islandsCount = 15,
                                             hierarchical = TRUE)
            
            oldTime <- Sys.time() 
            bbo_benchmark(optimizer, "hierarchy", name, budget = itersCount)
            newTime <- Sys.time() - oldTime
            print(newTime)
            cat("\n")
          }
        }
      }
    }
  }
}
