generate_data_SR <- function(n, gamma, model, treat_form, typeY, pi1, p_n_ratio){

  library(MASS)
  library(mvtnorm)

  alpha0 <- 0.15
  p0 <- ceiling(n * alpha0)
  beta0_full <- 1/(1:p0)^(1/4)*(-1)^c(1:p0)
  Sigma_true <- matrix(0, nrow = p0, ncol = p0)
  for (i in 1:p0) {
    for (j in 1:p0) {
      Sigma_true[i, j] <- 0.1 ** (abs(i - j))
    }
  }
  
  if(typeY != 'poisson'){
    X <- rmvt(n, sigma = Sigma_true, df = 3)
  }else{
    X0 <- rmvt(n, sigma = Sigma_true, df = 3)
    X <- pmin(pmax(X0, -3), 3)
    rm(X0)
  }
  
  beta <- gamma * beta0_full / norm(beta0_full,type='2')
  if(model == 'linear'){
    lp0 <- X %*% beta 
  }else if(model == 'nl1'){
    lp0 <- 1/2 * sign(X %*% beta) * (abs(X %*% beta))^(1/3) + cos(X %*% beta) + X %*% beta
  }
  
  if(treat_form == 'homo'){
    delta_X <- 1
  }else{
    delta_X <- 1 - 1/2 * pmin(X[, 1]^2, 5) + 1/4 * X[, 1:10] %*% beta[1:10]
  }
  lp1 <- lp0 + delta_X
  
  if (typeY == 'binomial') {
    r0 <- plogis(2 * lp0)
    r1 <- plogis(2 * lp1)
    Y1 <- rbinom(n, size = 1, prob = r1)
    Y0 <- rbinom(n, size = 1, prob = r0)

  }else if(typeY == 'poisson'){
    lp1_tran <- pmin(lp1, 4)
    lp0_tran <- pmin(lp0, 4)
    
    r1 <- exp(lp1_tran)
    r0 <- exp(lp0_tran)
    
    Y1 <- rpois(n,r1)
    Y0 <- rpois(n,r0)
  }else if(typeY == 'gaussian'){
    r1 <- lp1;
    r0 <- lp0
    Y1 <- r1 + rnorm(n)
    Y0 <- r0 + rnorm(n)
  }
  
  A <- rbinom(n, size = 1, prob = pi1)
  Y <- A * Y1 + (1 - A) * Y0
  p <- ceiling(round(n*p_n_ratio)) 
  if(p > ncol(X)){
    if(typeY != 'poisson'){
      X_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
    }else{
      X0_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
      X_noise <- pmin(pmax(X0_noise, -3), 3)
      rm(X0_noise)
    }
    X_obs <- cbind(X, X_noise)
  }else{
    X_obs <- X[, 1:p, drop = FALSE]
  }
  
  data_ls <- list(
    X = X_obs, Y= Y, A = A,
    Y1 = Y1, Y0 = Y0,
    r1 = r1, r0 = r0
  )
  return(data_ls)
}
