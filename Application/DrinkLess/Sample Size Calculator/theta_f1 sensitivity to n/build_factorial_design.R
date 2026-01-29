# This code is used to create factorial design with empty n
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  # m = c(30, 60),
                  m = c(30),
                  pt_shape = c("constant equal pk"),
                  # pt_shape = c("constant 0.5"),
                  ft_t_shape = c("constant", "linear_theta"),
                  ft_t_theta = seq(from = -1, to = 1, by = 0.01),
                  ft_t_theta2 = c(0),
                  gt_t_shape = c("constant", "quadratic_theta"),
                  gt_t_theta = c(0, 0.2), 
                  # for specifying gt when gt_t_shape is linear or quadratic
                  # ATE_t = c(1.2, 1.3, 1.4, 1.5),
                  ATE1_t = c(3.36),
                  ATE1_t_empirical = NA,
                  ATE2_t = c(3.6),
                  ATE2_t_empirical = NA,
                  ASPNC_t =c(0.08),
                  ASPNC_t_empirical = NA, 
                  #AvgTau_t = c(0.5, 0.8, 1),
                  AvgTau_t = c(1),
                  tau_t_shape = c("constant"),
                  beta1_t0 = NA,
                  beta1_t1 = NA,
                  beta1_t2 = NA,
                  beta2_t0 = NA,
                  beta2_t1 = NA,
                  beta2_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  ft_w_shape = NA,
                  ft_w_theta = NA,
                  ft_w_theta2 = NA,
                  gt_w_shape = NA,
                  gt_w_theta = 0,
                  beta1_w0 = NA,
                  beta1_w1 = NA,
                  beta1_w2 = NA,
                  beta2_w0 = NA,
                  beta2_w1 = NA,
                  beta2_w2 = NA,
                  se_beta1_w0 = NA,
                  se_beta2_w0 = NA,
                  alpha_w0 = NA,
                  alpha_w1 = NA,
                  alpha_w2 = NA,
                  ATE1_w = NA,
                  ATE1_w_empirical = NA,
                  ATE2_w = NA,
                  ATE2_w_empirical = NA,
                  ASPNC_w = NA,
                  ASPNC_w_empirical = NA,
                  AvgTau_w = NA,
                  tau_w_shape = NA,
                  p_w = NA,
                  q_w = NA,
                  n = NA,
                  p_star = NA,
                  error_indicator = 0,
                  stringsAsFactors = FALSE)

# remove simulation settings where pt*ft is not in the linear span of gt
SD <- SD[!(SD$ft_t_shape == "linear_theta" &
             SD$gt_t_shape == "constant"), ]
SD <- SD[!(SD$ft_t_shape == "quadratic_theta" &
             SD$gt_t_shape != "quadratic_theta"), ]

#remove constant with  theta not zero
SD <- SD[!(SD$ft_t_shape == "constant" &
             SD$ft_t_theta != 0), ]

# remove constant with theta2 not zero
SD <- SD[!(SD$ft_t_shape == "constant" &
             SD$ft_t_theta2 != 0), ]

SD <- SD[!(SD$gt_t_shape == "constant" &
             SD$gt_t_theta != 0), ]


SD <- SD[!(SD$AvgTau_t == 1 &
             SD$tau_t_shape != "constant"), ]

## add stuff that are known by simulation design
# taut
SD$AvgTau_w <- SD$AvgTau_t
SD$tau_w_shape <- SD$tau_t_shape

# gt
SD$gt_w_shape <- SD$gt_t_shape
SD$gt_w_theta <- SD$gt_t_theta
SD$ASPNC_w <- SD$ASPNC_t 

# ft
SD$ft_w_shape <- SD$ft_t_shape
SD$ft_w_theta <- SD$ft_t_theta
SD$ft_w_theta2 <- SD$ft_t_theta2
SD$ATE1_w <- SD$ATE1_t
SD$ATE2_w <- SD$ATE2_t

dim(SD)

saveRDS(SD, file = "factorial_design.RDS")
