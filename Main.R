# Calculate ATE with n=1e7-----------------------
N <- 1e7

library(MASS)
library(mvtnorm)
alpha0 <- 0.15
p0 <- ceiling(n * alpha0)
Sigma_true <- matrix(0,nrow=p0,ncol=p0)
for(i in 1:p0){
  for(j in 1:p0){
    Sigma_true[i,j] <- 0.1**(abs(i-j))
  }
}

set.seed(123)
beta0_full <- 1/(1:p0)^(1/4)*(-1)^c(1:p0)
ate_true_mat <- NULL
X <- rmvt(N, sigma = Sigma_true, df = 3)

for(typeY in c("binomial","gaussian","poisson")){
  if(typeY %in% c('poisson')){
    X <- pmin(pmax(X, -3), 3)
  }
  
  for(gamma in c(1)){ 
    beta <- gamma * beta0_full / norm(beta0_full,type='2')
    
    for(model in c('linear','nl1')){
      if(model == 'linear'){
        lp0 <- X %*% beta 
      }else if(model == 'nl1'){
        lp0 <- 1/2 * sign(X %*% beta) * (abs(X %*% beta))^(1/3) + cos(X %*% beta) + X %*% beta
      }
      for(treat_form in c('homo','hete')){
        if(treat_form == 'homo'){
          delta_X <- 1
        }else{
          delta_X <- 1 - 1/2 * pmin(X[, 1]^2, 5) + 1/4 * X[, 1:10] %*% beta[1:10]
        }
        lp1 <- lp0 + delta_X
        if (typeY == 'binomial') {
          r0 <- plogis(2 * lp0)
          r1 <- plogis(2 * lp1)
          Y1 <- rbinom(N, size = 1, prob = r1)
          Y0 <- rbinom(N, size = 1, prob = r0)
        }else if(typeY == 'poisson'){
          lp1_tran <- pmin(lp1, 4)
          lp0_tran <- pmin(lp0, 4)
          r1 <- exp(lp1_tran)
          r0 <- exp(lp0_tran)
          Y1 <- rpois(N,r1)
          Y0 <- rpois(N,r0)
        }else if(typeY == 'gaussian'){
          r1 <- lp1;
          r0 <- lp0
          Y1 <- r1 + rnorm(N)
          Y0 <- r0 + rnorm(N)
        }
        
        EY1 <- mean(Y1)
        EY0 <- mean(Y0)
        ate_true <- mean(Y1) - mean(Y0)
        ate_true_mat <- rbind(ate_true_mat,
                              c(gamma, model, treat_form, typeY, N, ate_true, EY0, EY1))
        
      }
    }
  }
}

colnames(ate_true_mat) <- c("gamma", "model", "treat_form", "typeY", "N", "ate_true", 'EY0', 'EY1')

ate_true_df <- as.data.frame(ate_true_mat)
ate_true_df$gamma <- as.numeric(ate_true_df$gamma)
ate_true_df$N <- as.numeric(ate_true_df$N)
ate_true_df$ate_true <- as.numeric(ate_true_df$ate_true)
ate_true_df$EY0 <- as.numeric(ate_true_df$EY0)
ate_true_df$EY1 <- as.numeric(ate_true_df$EY1)

require(dplyr)
ate_true_df <- ate_true_df %>%
  arrange(typeY) %>%
  as.data.frame()

saveRDS(ate_true_df, file = '../SimuData/ate_true_df.rds')



## robincar + gOB + cal  --------------------------

args_command <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args_command[1])
gamma <- as.numeric(args_command[2])
model <- args_command[3]
treat_form <- args_command[4]
typeY <- args_command[5]
p_n_ratio <- as.numeric(args_command[6])
rep_start <- as.numeric(args_command[7])

require(dplyr)
pi1_vec <- c(1/2, 1/3)
pi1_type <- 2
pi1 <- pi1_vec[pi1_type]
pi0 <- 1 - pi1

ate_true_df <- readRDS('../SimuData/ate_true_df.rds')
loc_data_type <- which(ate_true_df$gamma == gamma & 
                         ate_true_df$model == model &
                         ate_true_df$treat_form == treat_form &
                         ate_true_df$typeY == typeY)
ate_true <- ate_true_df$ate_true[loc_data_type]
EY1_true <- ate_true_df$EY1[loc_data_type]
EY0_true <- ate_true_df$EY0[loc_data_type]


times_mc_total <- 500
times_mc_each <- 50
rep_start_loc <- seq(rep_start,(rep_start + times_mc_total - 1),times_mc_each)

for(rep_start_sep in rep_start_loc){
  file_name <- paste0("result_robincar_n_", n,
                      '_alpha0_', 0.15,
                      "_gamma_", gamma,
                      '_model_', model,
                      '_treat_form_', treat_form,
                      "_typeY_", typeY,
                      "_p_n_ratio_", p_n_ratio,
                      '_pi1_type_', pi1_type,
                      '_rep_num(', rep_start_sep, '-', rep_start_sep + times_mc_each - 1, ').rds')
  if(!file.exists(paste0('../Results/robincar/',file_name))){
    require(doParallel); require(foreach)
    core_num <- 5
    type <- ifelse(.Platform$OS.type=='windows','PSOCK','FORK')
    cl <- makeCluster(core_num,type);  
    registerDoParallel(cl)
    set.seed(rep_start_sep)
    clusterSetRNGStream(cl, iseed = rep_start_sep)
    result.all.ls <- foreach(rep_num = rep_start_sep:(rep_start_sep + times_mc_each - 1))%dopar%{
      set.seed(rep_num)
      source('./fun_utils-simudata.R')
      data_ls <- generate_data_SR(n, gamma, model, treat_form, typeY, pi1, p_n_ratio)
      X <- data_ls$X; 
      A <- data_ls$A
      Y <- data_ls$Y
      
      n1 <- sum(A == 1)
      n0 <- sum(A == 0)
      p <- ncol(X)
      
      tau1_unadj <- mean(Y[A == 1]) 
      tau0_unadj <- mean(Y[A == 0])
      tau_unadj <- tau1_unadj - tau0_unadj
      
      var_tau1_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A))
      var_tau0_unadj <- 1 / length(A) * (var(Y[A==0]) / (1 - mean(A)))
      var_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A) + var(Y[A==0]) / (1 - mean(A)))
      
      
      source('./fun_base.R')
      result.base.ls <- fit.base(X, Y, A, typeY)
      result.base.ls
    }
    stopImplicitCluster();stopCluster(cl)
    
    
    tau_mat <- NULL
    var_tau_mat <- NULL
    
    for(i in 1:length(result.all.ls)){
      result.ls <- result.all.ls[[i]]
      tau_mat <- rbind(tau_mat,
                       c(true = ate_true, result.ls$point_vec))
      
      var_tau_mat <- rbind(var_tau_mat, result.ls$var_vec)
    }
    
    saveRDS(list(
      tau_mat = tau_mat,
      var_tau_mat = var_tau_mat
    ),file = paste0('../Results/robincar/', file_name))
  }
}



## adj2(-cal), adj2c(-cal)-----------------------------------
args_command <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args_command[1])
gamma <- as.numeric(args_command[2])
model <- args_command[3]
treat_form <- args_command[4]
typeY <- args_command[5]
p_n_ratio <- as.numeric(args_command[6])
rep_start <- as.numeric(args_command[7])

# n <- 400; gamma <- 1; model <- "linear"; treat_form <- 'homo'; typeY <- "gaussian"; p_n_ratio <- 0.05;
require(dplyr)
pi1_vec <- c(1/2, 1/3)
pi1_type <- 2
pi1 <- pi1_vec[pi1_type]
pi0 <- 1 - pi1


ate_true_df <- readRDS('../SimuData/ate_true_df.rds')
loc_data_type <- which(ate_true_df$gamma == gamma & 
                         ate_true_df$model == model &
                         ate_true_df$treat_form == treat_form &
                         ate_true_df$typeY == typeY)
ate_true <- ate_true_df$ate_true[loc_data_type]
EY1_true <- ate_true_df$EY1[loc_data_type]
EY0_true <- ate_true_df$EY0[loc_data_type]


times_mc_total <- 100
times_mc_each <- 50
rep_start_loc <- seq(rep_start, (rep_start + times_mc_total - 1), times_mc_each)

for(rep_start_sep in rep_start_loc){
  
  file_name <- paste0("result_adj2_adj2c_n_",n,
                      '_alpha0_', 0.15, 
                      "_gamma_", gamma, 
                      '_model_', model, 
                      '_treat_form_', treat_form, 
                      "_typeY_", typeY, 
                      "_p_n_ratio_", p_n_ratio, 
                      '_pi1_type_', pi1_type, 
                      '_rep_num(',rep_start_sep,'-',rep_start_sep + times_mc_each - 1,').rds')
  if(!file.exists(paste0('../Results/adj2_adj2c/', file_name))){
    require(doParallel); require(foreach)
    core_num <- 5
    type <- ifelse(.Platform$OS.type=='windows','PSOCK','FORK')
    cl <- makeCluster(core_num,type);  
    registerDoParallel(cl)
    set.seed(rep_start_sep)
    clusterSetRNGStream(cl, iseed = rep_start_sep)
    result.all.ls <- foreach(rep_num = rep_start_sep:(rep_start_sep + times_mc_each - 1))%dopar%{
      set.seed(rep_num) 
      source('./fun_utils-simudata.R')
      data_ls <- generate_data_SR(n, gamma, model, treat_form, typeY, pi1, p_n_ratio)
      X <- data_ls$X; 
      A <- data_ls$A
      Y <- data_ls$Y
      
      tau1_unadj <- mean(Y[A == 1]) 
      tau0_unadj <- mean(Y[A == 0])
      tau_unadj <- tau1_unadj - tau0_unadj
      
      var_tau1_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A))
      var_tau0_unadj <- 1 / length(A) * (var(Y[A==0]) / (1 - mean(A)))
      var_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A) + var(Y[A==0]) / (1 - mean(A)))
      
      result.unadj.ls <- list(
        tau = tau_unadj,
        var_tau = var_unadj,
        tau1 = tau1_unadj,
        var_tau1 = var_tau1_unadj,
        tau0 = tau0_unadj,
        var_tau0 = var_tau0_unadj
      )
      
      
      Xc <- cbind(1, scale(X, scale = FALSE))
      pi1_hat <- mean(A)
      
      require(HOIFCar)
      result.adj2.adj2c.ls <- fit.adj2.adj2c.Super(Y, Xc, A, intercept = TRUE, lc = FALSE)
      
      result.adj2.adj2c.ls_cal <- fit.adj2.adj2c.Super(Y, Xc, A, intercept = TRUE, lc = TRUE)
      result.adj2.adj2c.ls_cal <- lapply(result.adj2.adj2c.ls_cal, setNames, c('adj2_cal', 'adj2c_cal'))
      
      list(
        result.unadj.ls = result.unadj.ls,
        result.adj2.adj2c.ls = result.adj2.adj2c.ls,
        result.adj2.adj2c.ls_cal = result.adj2.adj2c.ls_cal
      )
    }
    stopImplicitCluster();stopCluster(cl)
    
    tau_mat <- NULL
    var_infl_tau_mat <- NULL
    for(i in 1:length(result.all.ls)){
      result.ls <- result.all.ls[[i]]
      result.unadj.ls <- result.ls$result.unadj.ls
      result.adj2.adj2c.ls <- result.ls$result.adj2.adj2c.ls
      result.adj2.adj2c.ls_cal <- result.ls$result.adj2.adj2c.ls_cal
      
      tau_mat <- rbind(tau_mat,
                       c(true = ate_true, 
                         unadj = result.unadj.ls$tau, 
                         result.adj2.adj2c.ls$tau_vec,
                         result.adj2.adj2c.ls_cal$tau_vec))
      
      var_infl_tau_mat <- rbind(var_infl_tau_mat, 
                              c(unadj = result.unadj.ls$var_tau, 
                                result.adj2.adj2c.ls$var_infl_vec,
                                result.adj2.adj2c.ls_cal$var_infl_vec))
    }
    
    saveRDS(list(
      tau_mat = tau_mat,
      var_infl_tau_mat = var_infl_tau_mat
    ),file = paste0('../Results/adj2_adj2c/', file_name))
    
  }
}


## JASA, JASA-cal---------------- 
args_command <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args_command[1])
gamma <- as.numeric(args_command[2])
model <- args_command[3]
treat_form <- args_command[4]
typeY <- args_command[5]
p_n_ratio <- as.numeric(args_command[6])
rep_start <- as.numeric(args_command[7])

# n <- 400; gamma <- 1; model <- "linear"; treat_form <- 'homo'; typeY <- "gaussian"; p_n_ratio <- 0.025; pi1_type <- 1
require(dplyr)
pi1_vec <- c(1/2, 1/3)
pi1_type <- 2
pi1 <- pi1_vec[pi1_type]
pi0 <- 1 - pi1

ate_true_df <- readRDS('../SimuData/ate_true_df.rds')
loc_data_type <- which(ate_true_df$gamma == gamma & 
                         ate_true_df$model == model &
                         ate_true_df$treat_form == treat_form &
                         ate_true_df$typeY == typeY)
ate_true <- ate_true_df$ate_true[loc_data_type]
EY1_true <- ate_true_df$EY1[loc_data_type]
EY0_true <- ate_true_df$EY0[loc_data_type]

times_mc_total <- 50
times_mc_each <- 50
rep_start_loc <- seq(rep_start,(rep_start + times_mc_total - 1),times_mc_each)

for(rep_start_sep in rep_start_loc){
  file_name <- paste0("result_jasa_n_", n,
                      '_alpha0_', 0.15,
                      "_gamma_", gamma,
                      '_model_', model,
                      '_treat_form_', treat_form,
                      "_typeY_", typeY,
                      "_p_n_ratio_", p_n_ratio,
                      '_pi1_type_', pi1_type,
                      '_rep_num(', rep_start_sep, '-', rep_start_sep + times_mc_each - 1, ').rds')
  if(!file.exists(paste0('../Results/jasa/', file_name))){
    require(doParallel); require(foreach)
    core_num <- 5
    type <- ifelse(.Platform$OS.type=='windows','PSOCK','FORK')
    cl <- makeCluster(core_num,type);  
    registerDoParallel(cl)
    set.seed(rep_start_sep)
    clusterSetRNGStream(cl, iseed = rep_start_sep)
    result.all.ls <- foreach(rep_num = rep_start_sep:(rep_start_sep + times_mc_each - 1))%dopar%{
      set.seed(rep_num) 
      source('./fun_utils-simudata.R')
      data_ls <- generate_data_SR(n, gamma, model, treat_form, typeY, pi1, p_n_ratio)
      X <- data_ls$X; 
      A <- data_ls$A
      Y <- data_ls$Y
      
      tau1_unadj <- mean(Y[A == 1]) 
      tau0_unadj <- mean(Y[A == 0])
      tau_unadj <- tau1_unadj - tau0_unadj
      
      var_tau1_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A))
      var_tau0_unadj <- 1 / length(A) * (var(Y[A==0]) / (1 - mean(A)))
      var_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A) + var(Y[A==0]) / (1 - mean(A)))
      
      result.unadj.ls <- list(
        tau = tau_unadj,
        var_tau = var_unadj,
        tau1 = tau1_unadj,
        var_tau1 = var_tau1_unadj,
        tau0 = tau0_unadj,
        var_tau0 = var_tau0_unadj
      )
      
      
      Xc <- cbind(1, scale(X, scale = FALSE))
      pi1_hat <- mean(A)
      
      require(HOIFCar)
      result.JASA_mu.ls <- fit.JASA(Y = Y, X = Xc, A = A, family = typeY, pi1 = pi1_hat,
                                    is.parallel = FALSE, opt_obj = "mu")
      result.JASA_mu.ls <- lapply(result.JASA_mu.ls, setNames, c('jasa(mu)', 'jasa(mu)_cal'))
      
      result.JASA_beta.ls <- fit.JASA(Y = Y, X = Xc, A = A, family = typeY, pi1 = pi1_hat,
                                    is.parallel = FALSE, opt_obj = "beta")
      result.JASA_beta.ls <- lapply(result.JASA_beta.ls, setNames, c('jasa(beta)', 'jasa(beta)_cal'))
      
      list(
        result.unadj.ls = result.unadj.ls,
        result.JASA_mu.ls = result.JASA_mu.ls,
        result.JASA_beta.ls = result.JASA_beta.ls
      )
    }
    stopImplicitCluster();stopCluster(cl)
    
    tau_mat <- NULL
    var_infl_tau_mat <- NULL
    result_cmb.JASA.ls <- vector('list', length = length(result.all.ls))
    for(i in 1:length(result.all.ls)){
      result.ls <- result.all.ls[[i]]
      result.unadj.ls <- result.ls$result.unadj.ls
      result.JASA_mu.ls <- result.ls$result.JASA_mu.ls
      result.JASA_beta.ls <- result.ls$result.JASA_beta.ls
      
      tau_mat <- rbind(tau_mat,
                       c(true = ate_true, unadj = result.unadj.ls$tau, 
                         result.JASA_mu.ls$tau_vec, result.JASA_beta.ls$tau_vec))
       
      var_infl_tau_mat <- rbind(var_infl_tau_mat, 
                                c(unadj = result.unadj.ls$var_tau, 
                                  result.JASA_mu.ls$var_infl_tau_vec, 
                                  result.JASA_beta.ls$var_infl_tau_vec))
      
      result_cmb.JASA.ls[[i]] <- list(
        result.JASA_mu.ls = result.ls$result.JASA_mu.ls,
        result.JASA_beta.ls = result.ls$result.JASA_beta.ls
      )
      
    }
    
    saveRDS(list(
      tau_mat = tau_mat,
      var_infl_tau_mat = var_infl_tau_mat,
      result_cmb.JASA.ls = result_cmb.JASA.ls
    ),file = paste0('../Results/jasa/', file_name))
    
  }
}


## LOO, LOO-Cal-In/Out ---------------------

args_command <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args_command[1])
gamma <- as.numeric(args_command[2])
model <- args_command[3]
treat_form <- args_command[4]
typeY <- args_command[5]
p_n_ratio <- as.numeric(args_command[6])
rep_start <- as.numeric(args_command[7])

# n <- 400; gamma <- 1; model <- "linear"; treat_form <- 'homo'; typeY <- "gaussian"; p_n_ratio <- 0.005; pi1_type <- 1
require(dplyr)
pi1_vec <- c(1/2, 1/3)
pi1_type <- 2
pi1 <- pi1_vec[pi1_type]
pi0 <- 1 - pi1

ate_true_df <- readRDS('../SimuData/ate_true_df.rds')
loc_data_type <- which(ate_true_df$gamma == gamma & 
                         ate_true_df$model == model &
                         ate_true_df$treat_form == treat_form &
                         ate_true_df$typeY == typeY)
ate_true <- ate_true_df$ate_true[loc_data_type]
EY1_true <- ate_true_df$EY1[loc_data_type]
EY0_true <- ate_true_df$EY0[loc_data_type]


times_mc_total <- 50
times_mc_each <- 50
rep_start_loc <- seq(rep_start,(rep_start + times_mc_total - 1),times_mc_each)

for(rep_start_sep in rep_start_loc){
  file_name <- paste0("result_loo_n_", n,
                      '_alpha0_', 0.15,
                      "_gamma_", gamma,
                      '_model_', model,
                      '_treat_form_', treat_form,
                      "_typeY_", typeY,
                      "_p_n_ratio_", p_n_ratio,
                      '_pi1_type_', pi1_type,
                      '_rep_num(', rep_start_sep, '-', rep_start_sep + times_mc_each - 1, ').rds')
  if(!file.exists(paste0('../Results/LOO/', file_name))){
    require(doParallel); require(foreach)
    core_num <- 5
    type <- ifelse(.Platform$OS.type=='windows', 'PSOCK', 'FORK')
    cl <- makeCluster(core_num, type);  
    registerDoParallel(cl)
    set.seed(rep_start_sep)
    clusterSetRNGStream(cl, iseed = rep_start_sep)
    result.all.ls <- foreach(rep_num = rep_start_sep:(rep_start_sep + times_mc_each - 1))%dopar%{
      set.seed(rep_num) 
      source('./fun_utils-simudata.R')
      data_ls <- generate_data_SR(n, gamma, model, treat_form, typeY, pi1, p_n_ratio)
      X <- data_ls$X; 
      A <- data_ls$A
      Y <- data_ls$Y
      
      n1 <- sum(A == 1)
      n0 <- sum(A == 0)
      p <- ncol(X)
      
      tau1_unadj <- mean(Y[A == 1]) 
      tau0_unadj <- mean(Y[A == 0])
      tau_unadj <- tau1_unadj - tau0_unadj
      
      var_tau1_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A))
      var_tau0_unadj <- 1 / length(A) * (var(Y[A==0]) / (1 - mean(A)))
      var_unadj <- 1 / length(A) * (var(Y[A==1]) / mean(A) + var(Y[A==0]) / (1 - mean(A)))
      
      result.unadj.ls <- list(
        tau = tau_unadj,
        var_tau = var_unadj,
        tau1 = tau1_unadj,
        var_tau1 = var_tau1_unadj,
        tau0 = tau0_unadj,
        var_tau0 = var_tau0_unadj
      )
      
      source('./fun_LeaveOneOut.R')
      pi1_hat <- mean(A)
      if(n1 > p & n0 > p){
        result.LOO.ls <- LeaveOneOut.Bern(observations = Y, covariates = X, allocation = A, 
                                          pi1 = pi1_hat, family = typeY)
      }else{
        tau1_vec <- tau0_vec <- tau_vec <- var_tau_vec <- rep(NA, 3)
        names(tau1_vec) <- names(tau0_vec) <- names(tau_vec) <- names(var_tau_vec) <- c('loo', 'loo_cal', 'loo_cal_out')
        y_hat_mat <- NULL
        result.LOO.ls <- list(
          tau0_vec = tau0_vec,
          tau1_vec = tau1_vec,
          tau_vec = tau_vec,
          var_tau_vec = var_tau_vec,
          y_hat_mat = y_hat_mat
        )
        
        
      }
      
      list(
        result.unadj.ls = result.unadj.ls,
        result.LOO.ls = result.LOO.ls
      )
    }
    stopImplicitCluster();stopCluster(cl)
    
    
    tau_mat <- NULL
    var_tau_mat <- NULL
    result_cmb.loo.ls <- vector('list', length = length(result.all.ls))
    for(i in 1:length(result.all.ls)){
      result.ls <- result.all.ls[[i]]
      result.unadj.ls <- result.ls$result.unadj.ls
      result.LOO.ls <- result.ls$result.LOO.ls
      
      tau_mat <- rbind(tau_mat,
                       c(true = ate_true, unadj = result.unadj.ls$tau, result.LOO.ls$tau_vec))
      
      var_tau_mat <- rbind(var_tau_mat, c(unadj = result.unadj.ls$var_tau, result.LOO.ls$var_tau_vec))
      
      result_cmb.loo.ls[[i]] <- result.ls$result.LOO.ls
      
    }
    
    
    saveRDS(list(
      tau_mat = tau_mat,
      var_tau_mat = var_tau_mat,
      result_cmb.LOO.ls = result_cmb.loo.ls
    ), file = paste0('../Results/LOO/', file_name))
    
  }
}



