# Code to calculate sample size 
# debugging from setting since it seems that the theoretical variance covariance is deifferent by scale of 5 
# this is the third version that allows different time pattern on availability. 
# nA is the number of Treatment option (control is included)
# p is the length of vector beta hat
# q is the length of vector alpha hat 
# nft is the length of moderator variable intercept is included
# pstar is the row length of contrast matrix, representing the rank of matrix L 
# f_t Defines marginal excursion effect MEE(t) under alternative together with beta. Assumed to be matrix of size t by nft.
# g_t Defines success probability null curve together with alpha. Assumed to be matrix of size T by q.(1 < nf < q)
# beta Length  p vector that defines marginal excursion effect MEE(t) under alternative together with f_t.
# alpha Length q vector that defines success probability null curve together with g_t.
# p_t Length T by nA vector of randomization probabilities at each time point
# L is contrast matrix with dimension pstar by p to compare marginal effect between treatment+
# sigma_matrix is the calculated variance matrix of c^t * beta
# gamma is type 1 error
# b is type II error
# tau is the vector that store average availability at time t E(I_t).
# if tau is non specified we assume always available at any point. 


compute_sigma <- function(f_t,
                          g_t,
                          alpha,
                          beta,            
                          p_t, 
                          tau
) {
  nA <- dim(pt)[2] # include ctr 
  p <- length(beta)
  ## The M and Sigma matrices (needed to compute lambda)
  M_matrix <- Sigma_matrix <- matrix(data=0, nrow=p, ncol=p) 

  t <- dim(f_t)[1]
  nft <- dim(f_t)[2]
  
  if(nft != p/(nA-1)){
    stop("Incorrect dimensions for f_t.")
  }
    
  # For each decision point 
  # (T = total # of decision points is taken as the length of p_t, for now)
  for (i in 1:t){
    
    # breaking down steps to identify bug and improve robustness of code
    this_f_t <- as.numeric(f_t[i, ])
    this_g_t <-  as.numeric(g_t[i, ])
    this_p_t <- as.numeric(p_t[i,])
    tau_t <- as.numeric(tau[i])
    

    # this_pt exclude assignment probability of no treatment
    this_pt <- this_p_t[2:nA]
    
    
    # calculate all the running sum constant for the matrices
    sum_ptfb <- 0  # calculate sum(pt*fbeta)
    sum_pe <- 0 # calculate e^(sum p_k e^-fbk -1)
    
    for (k in 1:length(this_pt)) {
      fb_k <- this_f_t %*% beta[((k-1) * length(this_f_t) + 1) :(k * length(this_f_t))]
      sum_ptfb <- sum_ptfb + (this_pt[k] * fb_k)[1]
      sum_pe <- sum_pe + this_pt[k] * (exp(-1* fb_k) -1)[1]
    }
    
    ut <- exp(sum_ptfb)
    
    # calculate SPNC 
    eY0 <-  exp(this_g_t %*% alpha)[1]
    
    # calculate f_t %*% f_t^T
    this_f_t_f_t <- outer(this_f_t, this_f_t)
    
    # M matrix 
    # create p_matrix
    p_matrix<- sigma_matrix <- matrix(data = 0, nrow = (nA-1), ncol = (nA-1))

    for (row in 1:(nA-1)) {
      for (col in 1:(nA-1)) {
        if (row == col) {
          # diagonal element of p_matrix
          p_matrix[row, col] = this_pt[row] * (1 - this_pt[row])
          # diagonal element of sigma matrix
          k <- row 
          fb_k <- (this_f_t %*% beta[((row-1) * length(this_f_t) + 1) :(row * length(this_f_t))])[1]
          # summation part
          sigma_matrix[row, col] = this_pt[row] * ( (1 - this_pt[row]) * (this_pt[row] + (1 - this_pt[row]) * exp(-1 * fb_k)  - eY0) +
                                                      this_pt[row] * (sum_pe - this_pt[row] * (exp(-1 * fb_k) - 1)) )
        }else{
          # off diagonal of p_matrix
          p_matrix[row, col] = -1 * this_pt[row] * this_pt[col]
          # off diagonal of sigma matrix
          k <- row
          j <- col
          fb_k <-  (this_f_t %*% beta[((k-1) * nft + 1) :(k * nft)])[1]
          fb_j <-  (this_f_t %*% beta[((j-1) * nft + 1) :(j * nft)])[1]
          # off diagonal of p_matrix
          sigma_matrix[row, col] =  this_pt[row] * this_pt[col] * ( 1- exp(-1 * fb_j)- exp(-1 * fb_k) + eY0 + sum_pe)
        }
      }
    } 

    # calculate the M matrix at time t 
    this_M <- kronecker(p_matrix, this_f_t_f_t) * tau_t * ut *eY0 # kroenecker product(see asymptotic notes)

    
    # calculate sigma matrix at time t 
    this_sigma <- kronecker(sigma_matrix, this_f_t_f_t) * tau_t * ut^2 *eY0 
    # running sum 
    M_matrix <- M_matrix + this_M
    Sigma_matrix <- Sigma_matrix + this_sigma
  }
  
  return(list(M_matrix = M_matrix, Sigma_matrix = Sigma_matrix))
  
}


compute_ncp <- function(x, beta, M_matrix,Sigma_matrix, L, tau){
# function to calculate the non centralize parameter
  if(det(M_matrix) == 0){
    stop("m_matrix must be nonsingular")
  }
  
  if(length(beta) != dim(M_matrix)[1]){
    stop("Dimensions of beta and m_matrix do not agree")
  }
  
  if(dim(M_matrix)[1] != dim(Sigma_matrix)[2] | 
     dim(M_matrix)[2] != dim(Sigma_matrix)[1]){
    stop("Dimensions of m_matrix and sigma_matrix do not agree")
  }
  
  if(det(solve(M_matrix) %*%
         Sigma_matrix %*%
         t(solve(M_matrix))) == 0){
    stop("sigma_matrix must be nonsingular")
  }
  
  return(as.numeric(x * t(L %*% beta) %*%
                      solve(L%*% solve(M_matrix)%*% Sigma_matrix %*% t(solve(M_matrix)) %*% t(L)) %*% 
                      L %*% beta))
}

mrt_mult_trt_ss <- function(f_t,
                            g_t, 
                            beta,
                            alpha,
                            p_t, 
                            gamma, 
                            b, 
                            L,
                            tau = rep(1,dim(f_t)[1]), 
                            exact = FALSE
){
  nA <- dim(pt)[2]
  t <- dim(f_t)[1]
  
  # degree of freedom of L*beta
  pstar <- dim(L)[1]
  # number of variable in g_t
  q <- dim(g_t)[2]

  varcov <- compute_sigma(f_t, g_t, alpha, beta, p_t, tau = tau)
  M_matrix <- varcov$M_matrix
  Sigma_matrix <- varcov$Sigma_matrix
  # Set up the function to solve using uniroot
  
  power_f <- function(n){
    
    right_hand_side <- pf(q=qf(p=(1-gamma), df1=pstar, df2=n-q-pstar), 
                          df1=pstar, 
                          df2=n-q-pstar, 
                          ncp=compute_ncp(n, beta,  M_matrix, Sigma_matrix, L, tau = tau))
    left_hand_side <- b
    return(right_hand_side - left_hand_side)
  }
  
  # find min n to achieve desired power
  sample_size <- uniroot(power_f, lower=pstar+q+1, upper=1000000)$root
  
  # round up if non-exact size is requested
  if(exact == FALSE){
    sample_size <- ceiling(sample_size)
  }
  
  return(sample_size)
}




#--------------------------------------------------------------------
# Test cases
#--------------------------------------------------------------------
# rm(list = ls(all = TRUE))
# source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/wcls_cat_trt_binary_outcome.R")
# source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/functions/utillity.R")
# 
# # two level treatment (ctrl, trt 1 , trt 2) with 15 time points
# m <- 15 # number of total Time
#prob_a <- c(0.2, 0.5, 0.3)
# 
# #simple case only intercept
# ## pt
# pt <- construct_pt(m, "constant equal pk")
# 
# ## taut_t
# taut_t <-
#   construct_taut_theta(m , 0.5, "constant")
# 
# ## gt_t, alpha_t
# gt_t <- construct_ftgt(m, "linear_theta")
# alpha_t <-
#   solve_alpha("linear_theta",
#               m,
#               ASPNC = 0.3,
#               gt_t,
#               taut_t,
#               theta = 0.3)
# 
# ## ft_t, beta1_t, beta2_t
# ft_t <- construct_ftgt(m, "constant")
# beta1_t <- solve_beta(
#   shape = "constant",
#   m = m,
#   ATE = 0.6,
#   ft = ft_t,
#   taut = taut_t
# )
# 
# beta2_t <- solve_beta(
#   "constant",
#   m,
#   ATE = 0.7,
#   ft = ft_t,
#   taut = taut_t
# )
# beta_t <- c(beta1_t, beta2_t)
# 
# L <- matrix(c(1,-1), nrow = 1)
# gamma <- 0.05
# b = 0.2
# 
# # build randomization probability matrix
# 
# mrt_mult_trt_ss(
#   f_t = ft_t,
#   g_t = gt_t,
#   beta = beta_t,
#   alpha = alpha_t,
#   p_t = pt,
#   gamma = gamma,
#   b = b,
#   L = L,
#   tau = taut_t
# )
# 
# 
# 
# mrt_mult_trt_ss(f_t, g_t, beta, alpha,  prob_a, mean_var, gamma, b, L, tau = tau)
# 
# simple case two covariates
## pt
# pt <- construct_pt(m, "constant equal pk")
# 
# ## taut_t
# taut_t <-
#   construct_taut_theta(m , 1, "constant")
# 
# ## gt_t, alpha_t
# gt_t <- construct_ftgt(m, "linear_theta")
# alpha_t <-
#   solve_alpha("linear_theta",
#               m,
#               ASPNC = 0.2,
#               gt_t,
#               taut_t,
#               theta = 0.3)
# 
# ## ft_t, beta1_t, beta2_t
# ft_t <- construct_ftgt(m, "linear_theta")
# beta1_t <- solve_beta(
#   shape = "linear_theta",
#   m = m,
#   ATE = 1.2,
#   ft = ft_t,
#   taut = taut_t,
#   theta = 0.3,
#   gt = gt_t,
#   alpha = alpha_t
# )
# beta2_t <- solve_beta(
#   "linear_theta",
#   m,
#   1.4,
#   ft = ft_t,
#   taut = taut_t,
#   theta = 0.3,
#   gt = gt_t,
#   alpha = alpha_t
# )
# beta_t <- c(beta1_t, beta2_t)
# 
# L <- matrix(cbind(diag(1, nrow = 2), -diag(1, nrow = 2)), nrow = 2)
# gamma <- 0.05
# b = 0.2
# 
# # build randomization probability matrix
# 
# mrt_mult_trt_ss(
#   f_t = ft_t,
#   g_t = gt_t,
#   beta = beta_t,
#   alpha = alpha_t,
#   p_t = pt,
#   gamma = gamma,
#   b = b,
#   L = L,
#   tau = taut_t
# )

# #--------------------------------------------------------------------
