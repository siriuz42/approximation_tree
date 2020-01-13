library(MASS)
source("calc_var.R")
ptestu <- function(test, model, correction = TRUE) {
  # Computes p-value for testing whether the means of all points in test are the same.
  #
  # Args:
  #   test: Covariates of test points
  #   model: Trained random forests model by library randomForest (set keep.inbag = TRUE)
  #   correction: Whether to apply bias correction for diagonal terms (variance)
  #
  # Returns:
  #   p-value 
  
  ans <- calc_var(test, model, correction)
  
  if (model$type == "regression") {
    
    y_pred <- ans$y_pred
    cov <- ans$var
    
    n_test <- dim(test)[1]
    
    df <- n_test
    
    I <- diag(n_test)
    J <- matrix(1, n_test, n_test)
    subu <- (I - J / n_test) %*% y_pred
    S <- (I - J / n_test) %*% cov %*% (I - J / n_test) 
    
    sta <- t(subu) %*% ginv(S) %*% subu
    
    return(pchisq(sta, df, lower.tail = FALSE))
  }
  
  if (model$type == "classification") {
    
    y_prob <- ans$y_prob
    cov <- ans$var
    
    n_test <- dim(test)[1]
    
    df <- n_test
    
    I <- diag(n_test)
    J <- matrix(1, n_test, n_test)
    subu <- (I - J / n_test) %*% y_prob
    S <- (I - J / n_test) %*% cov %*% (I - J / n_test) 
    
    sta <- t(subu) %*% ginv(S) %*% subu

    return(pchisq(sta, df, lower.tail = FALSE))
  }
  
}