library("foreach")
source("classicga.R")

islandsga <- function(fitness,
                      min,
                      max,
                      populationSize = 50,
                      islandsCount = 5,
                      migrationRate = 0.05,
                      migrationInterval = 10,
                      crossoverProb = 0.8,
                      mutationProb = 0.1,
                      elitismPercentage = 0.1,
                      iterations = 1000,
                      verbose = FALSE,
                      hierarchical = FALSE
) {
  # we consider only one specific hierarchical topology (in a form of binary tree)
  # 1 2 3 4 5 6 7 8
  #  9   A   B   C
  #    D       E
  #        F
  #
  if(hierarchical) {
    islandsCount <- 15
    hierarchyDirection <- c(9,9,10,10,11,11,12,12,13,13,14,14,15,15,NA)
  }
  
  suggestions <- matrix(nrow = 0, ncol = length(min))
  
  islSize <- max(1, floor(populationSize / islandsCount))
  migPop <- max(1, floor(migrationRate * islSize))
  numiter <- max(1, floor(iterations / migrationInterval))
  
  GAs <- vector(mode = "list", length = islandsCount)
  POPs <- rep(list(suggestions), times = islandsCount)
  for(iter in seq_len(numiter)) {
    GAs <- foreach(i. = seq_len(islandsCount)) %do% 
      { 
      classicga(
         fitness = fitness,
         min = min, 
         max = max,
         populationSize = islSize,
         crossoverProb = crossoverProb,
         mutationProb = mutationProb,
         iterations = migrationInterval,
         convergenceIters = migrationInterval,
         initialPopulation = POPs[[i.]]
         )
      }

    for(island in seq_len(islandsCount)) {
      # get summary of GAs evolution
      j <- seq((iter-1)*migrationInterval+1, iter*migrationInterval)
      
      # migration step
      if(!hierarchical) { 
        from <- island
        to <- (island %% islandsCount) + 1
      } else {
        from <- island
        to <- hierarchyDirection[[from]]
        if(verbose) print(paste0("Migrating from ", from, " to ", to))
        if(is.na(to)) {
          next
        }
      }
      
      ## select top individuals to migrationPop
      j <- order(GAs[[from]]$fitness, decreasing = TRUE)[seq(migPop)]
      migrationPop <- GAs[[from]]$population[j,,drop=FALSE]
      
      if(!hierarchical) {
        ## substitute random individuals
        j <- sample(islSize, size = migPop)
      } else {
        ## substitute random individuals but the elitist ones
        j <- sample(setdiff(seq(islSize),
                            order(GAs[[to]]$fitness, decreasing = TRUE)[seq(elitismPercentage)]),
                    size = migPop)
      }
      newpop <- rbind(GAs[[to]]$population[-j,,drop=FALSE], migrationPop)
      POPs[[to]] <- newpop
    }
  }
  
  fitnessValues <- lapply(GAs, function(x) max(x$fitness, na.rm = TRUE))
  # get islands' solutions
  solutions <- lapply(GAs, function(x) x$solution)
  fitnessValue <- max(unlist(fitnessValues), na.rm = TRUE)
  solution <- solutions[[which(fitnessValue == unlist(fitnessValues))[1]]]
  
  return(list(fitnessValues = fitnessValues, solutions = solutions, solution = solution))
}