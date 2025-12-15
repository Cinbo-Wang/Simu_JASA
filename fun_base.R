fit6.base <- function(X, Y, A, typeY, pi1 = NULL){
  require(RobinCar)
  n <- nrow(X); p <- ncol(X)
  n1 <- sum(A)
  n0 <- n-n1
  
  if(is.null(pi1)){
    pi1 <- mean(A)
  }
  pi0 <- 1 - pi1
  
  ### unadj: Different in means 
  tau_unadj <- mean(Y[A == 1]) - mean(Y[A == 0])
  var_unadj <- var(Y[A == 1]) / n1 + var(Y[A == 0]) / n0
  
  if(n1 < p | n0 < p){
    tau_anhecova <- NA
    tau_ob <- NA
    tau_cal <- NA
    tau_robincar_hete <- tau_robincar_hete_lc <- NA
    
    var_anhecova <- NA
    var_ob <- NA
    var_cal <- NA
    var_robincar_hete <- var_robincar_hete_lc <- NA
  }else{
    ### ANHECOVA: use gaussian 
    data_df <- data.frame(cbind(Y, A, X))
    colnames(data_df) <- c('Y', 'A', paste0('X', 1:ncol(X)))
    data_df$A <- as.factor(data_df$A)
    fit.anhecova <- robincar_linear(df = data_df,
                                    response_col = 'Y',
                                    treat_col = 'A',
                                    covariate_cols = paste0('X',1:ncol(X)),
                                    car_scheme = 'simple',
                                    adj_method = 'ANHECOVA',
                                    contrast_h = 'diff')
    tau_anhecova <- as.numeric(fit.anhecova$contrast$result$estimate)
    var_anhecova <- fit.anhecova$contrast$result$se^2
    
    ### Oaxaca-Blinder estimator 
    data_df <- data.frame(cbind(Y, X))
    colnames(data_df) <- c('Y', paste0('X', 1:ncol(X)))
    
    y1_hat <- tryCatch({
      fit1.glm <- glm(Y ~ .,
                      family = typeY,
                      data = data_df,
                      subset = (A == 1))
      y1_hat <- predict(fit1.glm, data_df, 'response')
      y1_hat
    }, error = function(e){
      rep(NA, length(Y))
    })
    y0_hat <- tryCatch({
      fit0.glm <- glm(Y ~ .,
                      family = typeY,
                      data = data_df,
                      subset = (A == 0))
      y0_hat <- predict(fit0.glm, data_df, 'response')
      y0_hat
    }, error = function(e){
      rep(NA, length(Y))
    })
   
    tau_ob <- mean(y1_hat - y0_hat)
    var_ob <- var((y1_hat - Y)[A == 1])/n1 + var((y0_hat - Y)[A == 0])/n0 + var(y1_hat - y0_hat)/n
    
   

    Xtilde <- data.frame(y1_hat,y0_hat)
    tau_cal_vec <- tryCatch({
      fit1.lm <- lm(Y~., data=Xtilde, subset = (A==1))
      fit0.lm <- lm(Y~., data=Xtilde, subset = (A==0))
      y1_tilde <- predict(fit1.lm, newdata = Xtilde)
      y0_tilde <- predict(fit0.lm, newdata = Xtilde)
      tau_cal <- mean(y1_tilde - y0_tilde)
      var_cal <- var((y1_tilde - Y)[A == 1])/n1 + var((y0_tilde - Y)[A == 0])/n0 + var(y1_tilde - y0_tilde)/n
      c(tau_cal, var_cal)
    }, error = function(e){
      c(NA, NA)
    })
    
    tau_cal <- tau_cal_vec[1]
    var_cal <- tau_cal_vec[2]
    
    ### Robincar: hete
    result_robincar_vec <- tryCatch({
      formula_hete <- paste0('Y ~ A*(', paste0('X', 1:ncol(X), collapse = " + "), ')')
      fit.hete <- robincar_glm(
        df = data_df,
        treat_col = 'A',
        response_col = 'Y',
        formula = formula_hete,
        car_scheme = 'simple',
        g_family = typeY
      )
      
      fit.hete.lc <- robincar_calibrate(result = fit.hete, joint=FALSE)
      fit.hete.diff <- robincar_contrast(fit.hete, contrast_h = 'diff')
      fit.hete.lc.diff <- robincar_contrast(fit.hete.lc, contrast_h = 'diff')
      
      # This is the same as Oaxaca-Blinder estimator
      tau_robincar_hete <- as.numeric(fit.hete.diff$result$estimate)
      # This is the same as No harm estimator
      tau_robincar_hete_lc <- as.numeric(fit.hete.lc.diff$result$estimate)
      
      # This is the same as Oaxaca-Blinder estimator
      var_robincar_hete <- as.numeric(fit.hete.diff$result$se^2)
      # This is the same as No harm estimator
      var_robincar_hete_lc <- as.numeric(fit.hete.lc.diff$result$se^2)
      
      c(tau_robincar_hete, tau_robincar_hete_lc, var_robincar_hete, var_robincar_hete_lc )
    }, error = function(e){
      c(NA, NA, NA, NA)
    })
    tau_robincar_hete <- result_robincar_vec[1]
    tau_robincar_hete_lc <- result_robincar_vec[2]
    var_robincar_hete <- result_robincar_vec[3]
    var_robincar_hete_lc <- result_robincar_vec[4]
  }
  
  
  point_base_vec <- c(tau_unadj,
                      tau_anhecova,
                      tau_robincar_hete,
                      tau_robincar_hete_lc,
                      tau_ob,
                      tau_cal)
  var_base_vec <- c(var_unadj,
                    var_anhecova,
                    var_robincar_hete,
                    var_robincar_hete_lc,
                    var_ob,
                    var_cal)
  names(point_base_vec) <- names(var_base_vec) <- c('unadj',
                                                    'ANHECOVA',
                                                    'robincar_hete',
                                                    'robincar_hete_lc',
                                                    'OB',
                                                    'cal')
  
  return(list(
    point_vec = point_base_vec,
    var_vec = var_base_vec
  ))
}
