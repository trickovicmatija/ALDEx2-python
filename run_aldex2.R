library(ALDEx2)


run_aldex <- function(counts, cond, mc.samples, test){
  #' Run all modules of ALDEx2 tool
  #'
  #'
  #' @param counts The dataframe containing counts with samples as columns 
  #' @param cond A list with len(num_of_samples) with conditions
  #' @param mc.samples Number of mc.samples
  #' @param test Type of statistical test
  #
  results <- aldex(counts, cond, mc.samples=mc.samples, test=test, effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
  
  return(results)
}

get_clr <- function(counts, cond, mc.samples,instance){
  #' Run clr module of ALDEx2 tool
  #'
  #'
  #' @param counts The dataframe containing counts with samples as columns 
  #' @param cond A list with len(num_of_samples) with conditions
  #' @param mc.samples Number of mc.samples
  #
  results_clr <- aldex.clr(counts, cond, mc.samples=mc.samples,, denom="all", verbose=FALSE)
  results <- getMonteCarloSample(results_clr,instance)
  return(results)
}