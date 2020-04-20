simulateSampling <- function(my.N, myprob) {
  vec <- rbinom(n = my.N, size = 1, prob = myprob)
  
  mdl <- glm(vec ~ 1, family = binomial)
  cis <- suppressMessages(plogis(confint(mdl)))
  
  out <- as.data.frame(t(cis))
  out$N <- my.N
  out$y <- plogis(coef(mdl))
  
  out
}

#' @param propN Number of simulations for same N in order to ascertain coverage.
#'              Defaults to 100.
simulateStudy <- function(N, myseed = NULL, myprob, propN = 100) {
  if (is.null(myseed)) {
    myseed <- sample(1:10^9, size = 1)
  }  
  set.seed(myseed)
  
  sims <- replicate(propN, simulateSampling(my.N = N, myprob = myprob),
                    simplify = FALSE)
  sims <- do.call(rbind, sims)
  colnames(sims) <- c("lci", "uci", "N", "y")
  
  sims$prob <- myprob
  sims$power <- ifelse(sims$prob >= sims$lci & sims$prob <= sims$uci,
                       yes = TRUE, no = FALSE)
  
  out <- colMeans(sims)
}

runBatch <- function(simseq, myprob, verbose = FALSE) {
  sims <- sapply(myprob, FUN = function(x) {
    sim <- sapply(simseq, FUN = simulateStudy, myprob = x, simplify = FALSE)
    sim <- do.call(rbind, sim)
    
    if (verbose) {
      print(sprintf("%s done", x))
    }
    
    sim
  }, simplify = FALSE)
  
  sims <- do.call(rbind, sims)
  sims <- as.data.frame(sims)
  
  sims$power_80 <- ifelse(sims$power >= 0.8, yes = TRUE, no = FALSE)
  sims$power_90 <- ifelse(sims$power >= 0.9, yes = TRUE, no = FALSE)
  sims$power_99 <- ifelse(sims$power >= 0.99, yes = TRUE, no = FALSE)
  
  sims
}