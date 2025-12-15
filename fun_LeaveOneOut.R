# Based on https://github.com/PeterLCohen/OaxacaBlinderCalibration/blob/main/Poisson_LeaveOneOut.R

LeaveOneOut.Bern = function(observations, covariates, allocation, 
                            pi1 = NULL, family = 'bernoulli'){  
  N <- length(allocation)
  nt <- sum(allocation)
  nc <- N - nt
  
  if(is.null(pi1)) pi1 <- mean(allocation)
  # U = (N/nt)*allocation - (N/nc)*(1 - allocation) # See https://arxiv.org/pdf/1708.01229.pdf equation 7 for definition
  # looVec = numeric(N)
  
  psi1CalVec  <- psi1Vec <- numeric(N)
  psi0CalVec  <- psi0Vec <- numeric(N)
  
  Y1_loo_hat <- Y0_loo_hat <- numeric(N)
  Y1_loo_cal_hat <- Y0_loo_cal_hat <- numeric(N)
  
  for(i in 1:N){
    # Hold out the ith data point
    Obs_holdOuti <- observations[-i]
    Cov_holdOuti <- covariates[-i,]
    treat_holdOuti <- allocation[-i]
    data <- data.frame(Y = Obs_holdOuti, X = Cov_holdOuti, treat = treat_holdOuti)
    
    # Fit on the data set without the ith data point
    mu1 <- glm(Y ~ . - treat, family = family, subset(data, treat == 1))
    mu0 <- glm(Y ~ . - treat, family = family, subset(data, treat == 0))
    # coef(mu1)
    
    tiHat = predict(mu1, data.frame(Y = observations[i],
                                    X = covariates[i, , drop = FALSE],
                                    treat = allocation[i]), "r")
    ciHat = predict(mu0, data.frame(Y = observations[i],
                                    X = covariates[i, , drop = FALSE],
                                    treat = allocation[i]), "r")
    Y1_loo_hat[i] <- tiHat
    Y0_loo_hat[i] <- ciHat
    # Gagnon-Barch & Wu's loo Estimator
    # mihat = (1 - (nt/N))*tiHat + (nt/N)*ciHat
    # looVec[i] = (observations[i] - mihat)*U[i]
    
    # The uncalibrated Poisson Rothe Leave-One-Out estimator
    # See the equation on the bottom of pg. 5 of http://www.christophrothe.net/papers/fca_apr2020.pdf
    # psiVec[i] = (tiHat - ciHat) +
    #   allocation[i]*(N/nt)*(observations[i] - tiHat) -
    #   (1 - allocation[i])*(N/nc)*(observations[i] - ciHat)
    
    psi1Vec[i] <- tiHat + allocation[i] / pi1 * (observations[i] - tiHat)
    psi0Vec[i] <- ciHat + (1 - allocation[i]) / (1 - pi1) * (observations[i] - ciHat)
    
    
    # Impute for all but the ith data point
    fittedT = predict(mu1, data, "r")
    fittedC = predict(mu0, data, "r")
    
    
    Xtilde = data.frame(fittedT, fittedC)
    XtildeT = Xtilde[treat_holdOuti==1,]
    XtildeC = Xtilde[treat_holdOuti==0,]
    Ytreat_minusi = Obs_holdOuti[treat_holdOuti==1]
    Ycontrol_minusi = Obs_holdOuti[treat_holdOuti==0]
    lmT2 = lm(Ytreat_minusi~., data = XtildeT)
    lmC2 = lm(Ycontrol_minusi~., data = XtildeC)
    
    ithPoint = data.frame(fittedT = tiHat, fittedC = ciHat)
    fittedT2 = predict(lmT2, newdata = ithPoint)
    fittedC2 = predict(lmC2, newdata = ithPoint)
    
    Y1_loo_cal_hat[i] <- fittedT2
    Y0_loo_cal_hat[i] <- fittedC2
    # See section 3.4 of http://www.christophrothe.net/papers/fca_apr2020.pdf
    # psiCalVec[i] = (fittedT2 - fittedC2) +
    #   allocation[i]*(N/nt)*(observations[i] - fittedT2) -
    #   (1 - allocation[i])*(N/nc)*(observations[i] - fittedC2)
    
    psi1CalVec[i] <- fittedT2 + allocation[i] / pi1 * (observations[i] - fittedT2)
    psi0CalVec[i] <- fittedC2 + (1 - allocation[i]) / (1 - pi1) * (observations[i] - fittedC2)
  }
  
  tau1_loo <- mean(psi1Vec) # Rothe Leave-One-Out estimator
  tau0_loo <- mean(psi0Vec)
  
  tau1_loo_cal <- mean(psi1CalVec)
  tau0_loo_cal <- mean(psi0CalVec)
  
  
  
  ### NEW: calibration out 
  Xtilde <- data.frame(y = observations, 
                       mu_c = Y0_loo_hat, 
                       mu_t = Y1_loo_hat)
  fit0 <- lm(y ~ ., data = Xtilde, subset = (allocation == 0))
  fit1 <- lm(y ~ ., data = Xtilde, subset = (allocation == 1))
  
  yc_hat_cal <- predict(fit0, newdata = Xtilde)
  yt_hat_cal <- predict(fit1, newdata = Xtilde)

  p <- ncol(covariates)
  if(p < 10){
    width <- max(observations) - min(observations)
    lb <- min(observations) - 0.1 * width; 
    ub <- max(observations) + 0.1 * width
  }else{
    lb <- -Inf; ub <- Inf
  }
  cond1 <- (min(yc_hat_cal) > lb & max(yc_hat_cal) < ub)
  cond2 <- (min(yt_hat_cal) > lb & max(yt_hat_cal) < ub)
  if(!(cond1 & cond2)){
    yc_hat_cal <- Y0_loo_hat
    yt_hat_cal <- Y1_loo_hat
  }
  
  tau1_loo_cal_out <- mean(allocation * observations / pi1 + (1 - allocation / pi1) * yt_hat_cal)
  tau0_loo_cal_out <- mean((1 - allocation) * observations / (1 - pi1) + (1 - (1 - allocation) / (1 - pi1)) * yc_hat_cal)
  
  tau_loo_cal_out <- tau1_loo_cal_out - tau0_loo_cal_out
  
  tau1_vec <- c(tau1_loo, tau1_loo_cal, tau1_loo_cal_out)
  tau0_vec <- c(tau0_loo, tau0_loo_cal, tau0_loo_cal_out)
  tau_vec <- tau1_vec - tau0_vec
  
  names(tau1_vec) <- names(tau0_vec) <- names(tau_vec) <- c('loo', 'loo_cal', 'loo_cal_out')
  
  
  est_MSE <- function(y,y_hat,t){
    diff_vec <- y - y_hat
    return(mean(diff_vec[allocation == t]^2))
  }
  
  Mt_hat <- est_MSE(observations, Y1_loo_hat, t=1)
  Mc_hat <- est_MSE(observations, Y0_loo_hat, t=0)
  var_loo <- 1 / N * ((1 - pi1) / pi1 * Mt_hat + pi1 / (1 - pi1) * Mc_hat + 2 * sqrt(Mt_hat * Mc_hat))
  
  
  Mt_hat <- est_MSE(observations, Y1_loo_cal_hat, t=1)
  Mc_hat <- est_MSE(observations, Y0_loo_cal_hat, t=0)
  var_loo_cal <- 1 / N * ((1 - pi1) / pi1 * Mt_hat + pi1 / (1 - pi1) * Mc_hat + 2 * sqrt(Mt_hat * Mc_hat))
  
  Mt_hat <- est_MSE(observations, yt_hat_cal, t=1)
  Mc_hat <- est_MSE(observations, yc_hat_cal, t=0)
  var_loo_cal_out <- 1 / N * ((1 - pi1) / pi1 * Mt_hat + pi1 / (1 - pi1) * Mc_hat + 2 * sqrt(Mt_hat * Mc_hat))
  
  
  var_tau_vec <- c(var_loo, var_loo_cal, var_loo_cal_out)
  names(var_tau_vec) <- c('loo', 'loo_cal', 'loo_cal_out')
  y_hat_mat <- cbind(Y0_loo_hat, Y1_loo_hat)
  return(list(
    tau0_vec = tau0_vec,
    tau1_vec = tau1_vec,
    tau_vec = tau_vec,
    var_tau_vec = var_tau_vec,
    y_hat_mat = y_hat_mat
  ))
  
}

