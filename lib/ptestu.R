library(MASS)

ptestu <- function(test, model) {
  
  ans = calc_var(test, model, correction = TRUE)
  
  if (model$type == "regression") {
    
    y_pred = ans$y_pred
    cov = ans$var
    
    n_test = dim(test)[1]
    
    df = n_test
    
    I <- diag(n_test)
    J <- matrix(1, n_test, n_test)
    subu <- (I - J / n_test) %*% y_pred
    S <- (I - J / n_test) %*% cov %*% (I - J / n_test) 
    
    sta <- t(subu) %*% ginv(S) %*% subu
    
    return(pchisq(sta, df, lower.tail = FALSE))
  }
  
  if (model$type == "classification") {
    
    y_prob = ans$y_prob
    cov = ans$var
    
    n_test = dim(test)[1]
    
    df = n_test
    
    I <- diag(n_test)
    J <- matrix(1, n_test, n_test)
    subu <- (I - J / n_test) %*% y_prob
    S <- (I - J / n_test) %*% cov %*% (I - J / n_test) 
    
    sta <- t(subu) %*% ginv(S) %*% subu

    return(pchisq(sta, df, lower.tail = FALSE))
  }
  
}