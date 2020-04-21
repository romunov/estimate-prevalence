#' From a simulated vector of positive and negative samples apply a
#' function which will simulate a rate of false negatives and positives.
#' 
#' I consider a false negative when subject is infected (1) but test
#' shows negative.
#' 
#' A false positive is a case when subject is not infected (0) but has
#' been diagnosed as such.
addFNFPbias <- function(x, fp, fn) {
  # Sorts values, zeros first and ones are clumped at the end.
  x <- sort(x)
  
  zeros <- which(x == 0)
  ones <- which(x == 1)
  
  # Work on false positives (0 -> 1).
  # If fp is small, almost no uninfected subjects will be tested as
  # false positive (basically no conversion from 0 to 1). Conversely,
  # if fp is large, almost all negative subjects will be tested as
  # false positive, in which case the test is a waste of resources.
  if (!is.na(fp)) {
    new.fp <- rbinom(n = length(zeros), size = 1, prob = fp)
    x[zeros] <- new.fp
  }
  
  # Work on false negative (1 -> 0).
  # By setting 1 - fn, we are removing the fn portion of ones from
  # otherwise positively tested pool. If fn is small, most cases that
  # have been tested positive will not be converted to false negative.
  # Vice versa, if fn is large, only a small fraction of positive cases
  # will remain positive and most will be converted to 0.
  if (!is.na(fn)) {
    new.fn <- rbinom(n = length(ones), size = 1, prob = 1 - fn)
    x[ones] <- new.fn
  }
  
  x
}

#' Run single simulation of sampling and estimation of coefficients.
simulateSampling <- function(my.N, myprob, fp, fn) {
  vec <- rbinom(n = my.N, size = 1, prob = myprob)
  
  vec <- addFNFPbias(x = vec, fp = fp, fn = fn)
  
  mdl <- glm(vec ~ 1, family = binomial)
  cis <- suppressMessages(plogis(confint(mdl)))
  
  out <- as.data.frame(t(cis))
  out$N <- my.N
  out$y <- plogis(coef(mdl))
  
  out
}

#' @param propN Number of simulations for same N in order to ascertain coverage.
#'              Defaults to 100.
simulateStudy <- function(N, myseed = NULL, myprob, propN, fp, fn) {
  if (is.null(myseed)) {
    myseed <- sample(1:10^9, size = 1)
  }  
  set.seed(myseed)
  
  sims <- replicate(propN, 
                    simulateSampling(my.N = N, myprob = myprob, fp = fp, 
                                     fn = fn),
                    simplify = FALSE)
  sims <- do.call(rbind, sims)
  colnames(sims) <- c("lci", "uci", "N", "y")
  
  sims$prob <- myprob
  sims$fp <- fp
  sims$fn <- fn
  sims$power <- ifelse(sims$prob >= sims$lci & sims$prob <= sims$uci,
                       yes = TRUE, no = FALSE)
  
  out <- colMeans(sims)
  out["sd"] <- sd(sims$y)
  
  out
}

runBatch <- function(simseq, myprob, propN = 100, fp = NA, fn = NA,
                     verbose = FALSE, cl) {
  clusterExport(cl = cl, varlist = c("simseq", "propN", "fp", "fn", "verbose"), 
                envir = environment())
  
  sims <- parSapply(cl = cl, X = myprob, FUN = function(x) {
    sim <- sapply(simseq, FUN = simulateStudy, 
                  myprob = x, 
                  propN = propN,
                  fp = fp,
                  fn = fn,
                  simplify = FALSE)
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

plotPower <- function(x) {
  ggplot(x, aes(x = N, y = power, color = as.factor(prob))) +
    theme_bw() +
    theme(legend.position = "top", legend.direction = "horizontal") +
    scale_color_brewer(palette = "Set1", name = "simulated prevalence") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(0, 10000, by = 2000)) +
    geom_vline(xintercept = 3000, color = "grey60", linetype = "dashed") +
    geom_line() +
    geom_point(size = 2.5)
}

plotCI <- function(x) {
  ggplot(x, aes(x = N, y = y)) +
    theme_bw() +
    theme(legend.position = "top") +
    geom_vline(xintercept = 3000, color = "grey60", linetype = "dashed") +
    geom_hline(aes(yintercept = prob), color = "grey60", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 10000, by = 2000)) +
    geom_pointrange(aes(ymin = y - sd, ymax = y + sd), size = 0.25) +
    facet_wrap(~ prob, scales = "free_y")
}