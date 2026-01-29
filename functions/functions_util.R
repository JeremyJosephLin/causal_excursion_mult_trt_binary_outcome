# These are codes needed to run for simulations 
wcls_summary <- function(fit){
  # create a function that calculate p-value of our wcls estimator
  # calculate degree of freedom 
  n <- fit$dims$n
  p <- fit$dims$p
  q <- fit$dims$q
  
  df = n - p - q 
  
  # extract beta 
  beta_hat <- fit$beta_hat
  
  # extract se betahat 
  se_beta <- fit$beta_se_adjusted
  
  # find critical value 
  t_val <- beta_hat / se_beta
  
  # pvalue is based on t distribution with n - p - q distribution 
  # find probability using t-distribution 
  pval <- 2*pt(abs(t_val), df= df, lower.tail = FALSE)
  
  coefficient <- data.frame(estimate = beta_hat , se = se_beta, t_val = t_val, pval = pval)
  colnames(coefficient) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  
  return(coefficient)
}


# between treatment marginal
wcls_glh <- function(fit, L){
  # q already include intercept
  q <- fit$dims$q
  n <- fit$dims$n
  p_star = dim(L)[1]
  
  Lbeta <- L %*%fit$beta_hat
  name_beta <- names(fit$beta_hat)
  
  # Variance covariance Matrix
  var_diff_adj <- L %*% fit$varcov_adjusted[name_beta, name_beta] %*% t(L)
  # standard error 
  se_Lbeta <-sqrt(diag(var_diff_adj)) 
  
  ### 1. calculate confidence interval
  # the distribution of Lbeta still will be asymptotically normal
  c <- qt(1 - 0.05/2, df = n - q - p_star)
  
  
  conf_int_adjusted <-
    cbind(Lbeta - c * se_Lbeta,
          Lbeta + c * se_Lbeta)
  
  colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")
  
  ### 2. Calculate p_value
  # find critical value 
  t_val <-  Lbeta /se_Lbeta
  
  # pvalue is based on t distribution with n - p - q distribution 
  # find probability using t-distribution 
  pval <- 2*pt(abs(t_val), df= n - q - p_star, lower.tail = FALSE)
  
  coefficient <- data.frame(estimate = Lbeta , se = se_Lbeta, t_val = t_val, pval = pval)
  colnames(coefficient) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  
  
  ### 3. Calculate p_value for simultaneous test 
  F_val <- (n- q- p_star + 1) / (p_star * (n - q)) * (t(Lbeta) %*% solve(var_diff_adj) %*% Lbeta)
  pval_simultaneous <- pf(F_val, df1 = p_star, df2 = n - p_star - q + 1, lower.tail = FALSE)
  
  
  result <- list(conf_int = conf_int_adjusted, summary = coefficient, simultaneous_test = pval_simultaneous)
  
  return(result)
}
