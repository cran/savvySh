## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(knitr)

# Parameters
set.seed(123)
n_val <- 1000
p_val <- 10
rho_val <- 0.75
sigma_val <- 5
mu_val <- 0

# Correlation matrix
sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

# True beta
theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

# Simulate data
Sigma <- sigma.rho(rho_val, p_val)
X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
X_intercept <- cbind(1, X)
beta_true <- theta_func(p_val + 1)
y <- rnorm(n_val, mean = X_intercept %*% beta_true, sd = sigma_val)

# Fit models
ols_fit <- lm(y ~ X)
beta_ols <- coef(ols_fit)

multi_results <- savvySh(X, y, model_class = "Multiplicative", include_Sh = TRUE)
beta_St  <- coef(multi_results, "St")
beta_DSh <- coef(multi_results, "DSh")
beta_Sh  <- coef(multi_results, "Sh")

slab_results <- savvySh(X, y, model_class = "Slab")
beta_SR  <- coef(slab_results, "SR")
beta_GSR <- coef(slab_results, "GSR")

## -----------------------------------------------------------------------------
# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "St", "DSh", "Sh", "SR", "GSR"),
  L2_Distance = c(
    sqrt(sum((beta_ols - beta_true)^2)),
    sqrt(sum((beta_St  - beta_true)^2)),
    sqrt(sum((beta_DSh - beta_true)^2)),
    sqrt(sum((beta_Sh  - beta_true)^2)),
    sqrt(sum((beta_SR  - beta_true)^2)),
    sqrt(sum((beta_GSR - beta_true)^2))
  )
)

kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_ols, 3),
  St = round(beta_St, 3),
  DSh = round(beta_DSh, 3),
  Sh = round(beta_Sh, 3),
  SR = round(beta_SR, 3),
  GSR = round(beta_GSR, 3),
  True = round(beta_true, 3)
)

kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(knitr)

# Parameters
set.seed(123)
n_val <- 1000
p_val <- 10
rho_val <- 0.75
df_val = 50/24 
mu_val <- 0

# Correlation matrix
sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

# True beta
theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

# Simulate data
Sigma <- sigma.rho(rho_val, p_val)
X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
X_intercept <- cbind(1, X)
beta_true <- theta_func(p_val + 1)
y <- as.vector(X_intercept %*% beta_true) + rt(n = n_val, df = df_val)

# Fit models
ols_fit <- lm(y ~ X)
beta_ols <- coef(ols_fit)

multi_results <- savvySh(X, y, model_class = "Multiplicative", include_Sh = TRUE)
beta_St  <- coef(multi_results, "St")
beta_DSh <- coef(multi_results, "DSh")
beta_Sh  <- coef(multi_results, "Sh")

slab_results <- savvySh(X, y, model_class = "Slab")
beta_SR  <- coef(slab_results, "SR")
beta_GSR <- coef(slab_results, "GSR")

## -----------------------------------------------------------------------------
# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "St", "DSh", "Sh", "SR", "GSR"),
  L2_Distance = c(
    sqrt(sum((beta_ols - beta_true)^2)),
    sqrt(sum((beta_St  - beta_true)^2)),
    sqrt(sum((beta_DSh - beta_true)^2)),
    sqrt(sum((beta_Sh  - beta_true)^2)),
    sqrt(sum((beta_SR  - beta_true)^2)),
    sqrt(sum((beta_GSR - beta_true)^2))
  )
)

kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_ols, 3),
  St = round(beta_St, 3),
  DSh = round(beta_DSh, 3),
  Sh = round(beta_Sh, 3),
  SR = round(beta_SR, 3),
  GSR = round(beta_GSR, 3),
  True = round(beta_true, 3)
)

kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(knitr)

# Parameters
set.seed(123)
n_val <- 1000
p_val <- 10
rho_val <- 0.75
q_0 <- 0.01
sigma_val <- 5

# Correlation matrix
sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

sigma.temp <- sigma.rho(rho_val, p_val) 
Z <- mvrnorm(n_val, mu = rep(0, p_val), Sigma = sigma.temp) 
X <- apply(Z, 2, function(z_col) qbinom(pnorm(z_col), size = 2, prob = q_0))  
X_intercept <- cbind(1, X)
beta_true <- c(0, runif(p_val, 0.01, 0.3))
y <- rnorm(n_val, mean = as.vector(X_intercept %*% beta_true), sd = sigma_val)

# Fit models
ols_fit <- lm(y ~ X)
beta_ols <- coef(ols_fit)

multi_results <- savvySh(X, y, model_class = "Multiplicative", include_Sh = TRUE)
beta_St  <- coef(multi_results, "St")
beta_DSh <- coef(multi_results, "DSh")
beta_Sh  <- coef(multi_results, "Sh")

slab_results <- savvySh(X, y, model_class = "Slab")
beta_SR  <- coef(slab_results, "SR")
beta_GSR <- coef(slab_results, "GSR")

## -----------------------------------------------------------------------------
# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "St", "DSh", "Sh", "SR", "GSR"),
  L2_Distance = c(
    sqrt(sum((beta_ols - beta_true)^2)),
    sqrt(sum((beta_St  - beta_true)^2)),
    sqrt(sum((beta_DSh - beta_true)^2)),
    sqrt(sum((beta_Sh  - beta_true)^2)),
    sqrt(sum((beta_SR  - beta_true)^2)),
    sqrt(sum((beta_GSR - beta_true)^2))
  )
)

kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_ols, 3),
  St = round(beta_St, 3),
  DSh = round(beta_DSh, 3),
  Sh = round(beta_Sh, 3),
  SR = round(beta_SR, 3),
  GSR = round(beta_GSR, 3),
  True = round(beta_true, 3)
)

kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(knitr)

# Parameters
set.seed(1)
n_val <- 1000
p_val <- 10
rho_val <- 0.75
sigma_val <- 5
mu_val <- 0

# Correlation matrix
sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

# True beta
theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

# Simulate data
Sigma <- sigma.rho(rho_val, p_val)
X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
X_centred <- scale(X, center = TRUE, scale = FALSE)
beta_true <- theta_func(p_val)
y <- rnorm(n_val, mean = X_centred %*% beta_true, sd = sigma_val)
y_centred <- scale(y, center = TRUE, scale = FALSE)

# Fit models
ols_fit <- lm(y_centred ~ X_centred-1)
beta_ols <- coef(ols_fit)

linear_results <- savvySh(X_centred, y_centred, model_class = "Linear")
beta_LSh <- coef(linear_results, "LSh")

## -----------------------------------------------------------------------------

# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "LSh"),
  L2_Distance = c(
    sqrt(sum((beta_ols - beta_true)^2)),
    sqrt(sum((beta_LSh  - beta_true)^2))
  )
)

kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_ols, 3),
  LSh = round(beta_LSh, 3),
  True = round(beta_true, 3)
)

kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(knitr)

# Parameters
set.seed(123)
n_val <- 1000
p_val <- 10
df_val = 50/24 
sigma_val <- 5
mu_val <- 0

# Correlation matrix
sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

# True beta
theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

# Simulate data
Sigma <- sigma.rho(rho_val, p_val)
X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
X_centred <- scale(X, center = TRUE, scale = FALSE)
beta_true <- theta_func(p_val)
y <- as.vector(X_centred %*% beta_true) + rt(n = n_val, df = df_val)
y_centred <- scale(y, center = TRUE, scale = FALSE)

# Fit models
ols_fit <- lm(y_centred ~ X_centred-1)
beta_ols <- coef(ols_fit)

linear_results <- savvySh(X_centred, y_centred, model_class = "Linear")
beta_LSh <- coef(linear_results, "LSh")

## -----------------------------------------------------------------------------
# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "LSh"),
  L2_Distance = c(
    sqrt(sum((beta_ols - beta_true)^2)),
    sqrt(sum((beta_LSh  - beta_true)^2))
  )
)

kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_ols, 3),
  LSh = round(beta_LSh, 3),
  True = round(beta_true, 3)
)

kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## -----------------------------------------------------------------------------
# Load packages
library(savvySh)
library(MASS)
library(glmnet)
library(knitr)

# Parameters
set.seed(123)
n_val <- 1000
p_val <- 10
f_val <- 5
sigma_val <- 5

# True beta
theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

A <- matrix(rnorm(n_val * f_val), nrow = n_val)  
Z <- matrix(rnorm(f_val * p_val), nrow = f_val) 
X <- A %*% Z  # n x p matrix
noise <- matrix(rnorm(n_val * p_val, sd = sqrt(10^(-6))), nrow = n_val)
X_noisy <- X + noise
X_intercept <- cbind(rep(1, n_val), X_noisy)
beta_true <- theta_func(p_val + 1)
y <- rnorm(n_val,mean=as.vector(X_intercept%*%beta_true),sd=sigma_val)

# Fit models
OLS_results <- lm(y~X_noisy)
beta_OLS <- coef(OLS_results)

glmnet_fit <- cv.glmnet(X_noisy, y, alpha = 0)
lambda_min_RR_glmnet <- glmnet_fit$lambda.min
beta_RR <- as.vector(coef(glmnet_fit, s = "lambda.min"))

SRR_results <- savvySh(X_noisy, y, model_class = "ShrinkageRR")
beta_SRR <- coef(SRR_results, "SRR")

## -----------------------------------------------------------------------------
# L2 comparison
l2_table <- data.frame(
  Method = c("OLS", "RR", "SRR"),
  L2_Distance = c(
     sqrt(sum((beta_OLS - beta_true)^2)),
    sqrt(sum((beta_RR - beta_true)^2)),
    sqrt(sum((beta_SRR  - beta_true)^2))
  )
)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients")

# Coefficient comparison table
coef_table <- data.frame(
  OLS = round(beta_OLS, 3),
  RR = round(beta_RR, 3),
  SRR = round(beta_SRR, 3),
  True = round(beta_true, 3)
)
kable(coef_table, caption = "Estimated Coefficients by Method (rounded)")

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("Ziwei-ChenChen/savvyGLM")
#  library(savvyGLM)
#  library(savvySh)
#  library(savvyGLM)
#  library(MASS)
#  
#  standardize_features<-function(dataset){
#    dataset[2:(ncol(dataset)-1)] <- as.data.frame(scale(dataset[2:(ncol(dataset)-1)]))
#    return(dataset)
#  }
#  
#  set_classes<-function(data){
#    data[,ncol(data)]<-replace(data[,ncol(data)], data[,ncol(data)]<1, 0)
#    data[,ncol(data)]<-replace(data[,ncol(data)], data[,ncol(data)] %in% c(1,2,2.5,3,3.5), 1)
#    data[,ncol(data)]<-replace(data[,ncol(data)], data[,ncol(data)] %in% c(4,5,6), 2)
#    data[,ncol(data)]<-replace(data[,ncol(data)], data[,ncol(data)]>=7, 3)
#    return(data)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  fit_and_return_coefficients <- function(x, y) {
#  
#    control_list <- list(maxit = 200, epsilon = 1e-6, trace = TRUE)
#    family_type <- poisson(link = "log")
#  
#    # Fitting models
#    opt_glm2_OLS <- glm.fit2(x, y, family = family_type, control = control_list)
#    opt_glm2_SR <- savvy_glm.fit2(x, y, model_class = "SR", family = family_type, control = control_list)
#    opt_glm2_GSR <- savvy_glm.fit2(x, y, model_class = "GSR", family = family_type, control = control_list)
#    opt_glm2_St <- savvy_glm.fit2(x, y, model_class = "St", family = family_type, control = control_list)
#    opt_glm2_DSh <- savvy_glm.fit2(x, y, model_class = "DSh", family = family_type, control = control_list)
#    opt_glm2_Sh <- savvy_glm.fit2(x, y, model_class = "Sh", family = family_type, control = control_list)
#  
#    return(list(
#      glm2_OLS_result = opt_glm2_OLS$coefficients,
#      glm2_SR_result = opt_glm2_SR$coefficients,
#      glm2_GSR_result = opt_glm2_GSR$coefficients,
#      glm2_St_result = opt_glm2_St$coefficients,
#      glm2_DSh_result = opt_glm2_DSh$coefficients,
#      glm2_Sh_result = opt_glm2_Sh$coefficients
#    ))
#  }
#  
#  test_model <- function(glm_coefficients, data_X, data_Y) {
#    upper_limit <- 3  # =3 for 4 classes; =9 for 10 classes
#  
#    ### Model 1 ---> OLS ###
#    predicted_glm2_OLS <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_OLS_result)))
#    predicted_glm2_OLS <- ifelse(predicted_glm2_OLS <= upper_limit, predicted_glm2_OLS, upper_limit)
#  
#    ### Model 2 ---> SR ###
#    predicted_glm2_SR <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_SR_result)))
#    predicted_glm2_SR <- ifelse(predicted_glm2_SR <= upper_limit, predicted_glm2_SR, upper_limit)
#  
#    ### Model 3 ---> GSR ###
#    predicted_glm2_GSR <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_GSR_result)))
#    predicted_glm2_GSR <- ifelse(predicted_glm2_GSR <= upper_limit, predicted_glm2_GSR, upper_limit)
#  
#    ### Model 4 ---> Stein ###
#    predicted_glm2_St <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_St_result)))
#    predicted_glm2_St <- ifelse(predicted_glm2_St <= upper_limit, predicted_glm2_St, upper_limit)
#  
#    ### Model 5 ---> Diagonal Shrinkage ###
#    predicted_glm2_DSh <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_DSh_result)))
#    predicted_glm2_DSh <- ifelse(predicted_glm2_DSh <= upper_limit, predicted_glm2_DSh, upper_limit)
#  
#    ### Model 6 ---> Sh ###
#    predicted_glm2_Sh <- floor(exp(data_X %*% as.matrix(glm_coefficients$glm2_Sh_result)))
#    predicted_glm2_Sh <- ifelse(predicted_glm2_Sh <= upper_limit, predicted_glm2_Sh, upper_limit)
#  
#    print(max(predicted_glm2_OLS))
#    print(max(predicted_glm2_SR))
#  
#    r_OLS <- c(mse(data_Y, predicted_glm2_OLS), rmse(data_Y, predicted_glm2_OLS), mae(data_Y, predicted_glm2_OLS))
#    r_SR <- c(mse(data_Y, predicted_glm2_SR), rmse(data_Y, predicted_glm2_SR), mae(data_Y, predicted_glm2_SR))
#    r_GSR <- c(mse(data_Y, predicted_glm2_GSR), rmse(data_Y, predicted_glm2_GSR), mae(data_Y, predicted_glm2_GSR))
#    r_St <- c(mse(data_Y, predicted_glm2_St), rmse(data_Y, predicted_glm2_St), mae(data_Y, predicted_glm2_St))
#    r_DSh <- c(mse(data_Y, predicted_glm2_DSh), rmse(data_Y, predicted_glm2_DSh), mae(data_Y, predicted_glm2_DSh))
#    r_Sh <- c(mse(data_Y, predicted_glm2_Sh), rmse(data_Y, predicted_glm2_Sh), mae(data_Y, predicted_glm2_Sh))
#  
#    return(list(
#      results_OLS = r_OLS,
#      results_SR = r_SR,
#      results_GSR = r_GSR,
#      results_St = r_St,
#      results_DSh = r_DSh,
#      results_Sh = r_Sh
#    ))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  out_of_sample_performance <- function(percentage = 0.3, no_trials = 50, filein = "data_for_regression10.csv", fileout = "results10.csv") {
#  
#    data <- read.csv(filein)
#    agregated_results<-rep(0,4*3)
#  
#    data<-set_classes(data)   # for 4 classes
#    #data[,ncol(data)]<-floor(data[,ncol(data)])  # for 10 classes -- the scores 2.5 and 3.5 are floored
#  
#    data<-standardize_features(data)
#  
#    for (r in 1:no_trials){
#      #for train-test split with percentage e.g. 70-30
#      bound <- floor(nrow(data)*(1-percentage)) #define % of training
#      data <- data[sample(nrow(data)), ]           #sample rows
#      train <- data[1:bound, ]
#      test <- data[(bound+1):nrow(data), ]
#  
#      #training
#      X.tilde <- train[,-ncol(train)]
#      X <-  X.tilde[,-1]
#      train_X <- as.matrix(train[,-ncol(train)])
#      train_Y <-  train[,ncol(train)]
#      glm_coefficients<-fit_and_return_coefficients(train_X, train_Y)
#  
#      #test
#      test_X <- as.matrix(test[,-ncol(test)])
#      test_Y <-  as.matrix(test[,ncol(test)])
#      results<-test_model(glm_coefficients, test_X, test_Y)
#      agregated_results<-agregated_results+unlist(results)
#    }
#  
#    # Average results over trials
#    aggregated_results <- aggregated_results / no_trials
#    df <- data.frame(matrix(unlist(agregated_results), nrow=length(results), byrow=TRUE))
#    colnames(df) <- c("MSE", "RMSE", "MAE")
#    rownames(df) <- c("GLM2", "SR", "GSR", "St", "DSh", "Sh")
#    write.csv(df, fileout)
#  }
#  
#  input_file_path <- "data_for_regression10.csv"
#  output_file_path <- "results10.csv"
#  run_performance_test <- out_of_sample_performance(
#    percentage = 0.3,
#    no_trials = 50,
#    filein = input_file_path,
#    fileout = output_file_path
#  )

## ----eval=FALSE---------------------------------------------------------------
#  library(savvySh)
#  library(MASS)
#  library(glmnet)
#  library(PerformanceAnalytics)
#  library(lubridate)
#  library(quadprog)
#  library(xts)
#  library(POET)
#  
#  data <- returns_441
#  data$Date <- as.Date(as.character(data$Date), format = "%Y%m%d")
#  colnames(data)[2:442] <- paste0("Company", 1:441)
#  
#  training_size <- 5 * 252
#  testing_size  <- 3 * 21
#  step_size <- 3 * 21
#  n_total <- nrow(data)
#  max_windows <- floor((n_total - training_size - testing_size) / step_size) + 1
#  cat("Total rows:", n_total, "\n")
#  cat("Max windows:", max_windows, "\n")
#  
#  get_full_weights <- function(est_vector) {
#    w <- est_vector[-1]
#    w_last <- 1 - sum(w)
#    return(c(w, w_last))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  poet_est_99 <- function(x, y) {
#    x <- as.matrix(x)
#    y <- as.numeric(y)
#    n <- nrow(x)
#    p <- ncol(x)
#  
#    x_mean <- colMeans(x)
#    y_mean <- mean(y)
#    x_c <- x - matrix(rep(x_mean, each = n), n, p)
#    y_c <- y - y_mean
#    Y=t(x_c)
#  
#    # Choose K to explain 99% variance
#    eigvals <- eigen(cov(x_c), symmetric = TRUE, only.values = TRUE)$values
#    eigvals_sorted <- sort(eigvals, decreasing = TRUE)
#    cumsum_vals <- cumsum(eigvals_sorted) / sum(eigvals_sorted)
#    K <- which(cumsum_vals >= 0.99)[1]
#    poet_result <- POET::POET(Y, K = K)
#    Sigma_hat <- poet_result$SigmaY
#  
#    xcyc <- crossprod(x_c, y_c) / n
#    beta_1p <- as.numeric(solve(Sigma_hat, xcyc))
#    beta_0 <- y_mean - sum(x_mean * beta_1p)
#    beta_full <- c(beta_0, beta_1p)
#  
#    return(as.vector(beta_full))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  rolling_annual_expected_returns <- data.frame()
#  rolling_annual_sharpe_ratios <- data.frame()
#  rolling_annual_volatilities <- data.frame()
#  
#  for (window_index in seq_len(max_windows)) {
#    start_index <- 1 + (window_index - 1) * step_size
#    train_start <- start_index
#    train_end <- start_index + training_size - 1
#    test_start <- train_end + 1
#    test_end <- train_end + testing_size
#  
#    train_data <- data[train_start:train_end, ]
#    test_data <- data[test_start:test_end, ]
#    train_returns <- as.matrix(train_data[, -1])
#    test_returns <- as.matrix(test_data[, -1])
#  
#    # Create X_train and Y_train
#    # Y_train is last column (Company441)
#    # X_train is (Y_train - first 440 columns)
#    Y_train <- train_returns[, 441]
#    X_train <- matrix(Y_train, nrow = nrow(train_returns), ncol = 440) - train_returns[, 1:440]
#  
#    # Fit all estimators
#    est_results <- list(
#      OLS = OLS_est(X_train, Y_train)$est,
#      RR = RR_est(X_train, Y_train)$est,
#      POET_99 = poet_est_99(X_train, Y_train),
#      St = St_ost(X_train, Y_train),
#      DSh = DSh_ost(X_train, Y_train),
#      Sh = Sh_ost(X_train, Y_train),
#      SR = SR_ost(X_train, Y_train),
#      GSR = GSR_ost(X_train, Y_train),
#      SRR = SRR_shrinkage_ost(X_train, Y_train)
#    )
#    weights_list <- lapply(est_results, get_full_weights)
#    names(weights_list) <- names(est_results)
#  
#    test_dates <- as.Date(test_data$Date)
#    test_returns_xts <- xts(test_returns, order.by = test_dates)
#    daily_returns_list <- lapply(weights_list, function(w) {
#      rp <- Return.portfolio(R = test_returns_xts, weights = w)
#      return(as.numeric(rp))
#    })
#    daily_values_list <- lapply(daily_returns_list, function(r) {
#      cum_val <- cumprod(1 + r)
#      return(cum_val)
#    })
#  
#    model_names <- names(daily_returns_list)
#    n_test <- length(test_start:test_end)
#    daily_values_mat <- matrix(0, nrow = length(model_names), ncol = n_test)
#    daily_returns_mat <- matrix(0, nrow = length(model_names), ncol = n_test)
#    rownames(daily_values_mat) <- model_names
#    rownames(daily_returns_mat) <- model_names
#    for (i in seq_along(model_names)) {
#      daily_values_mat[i, ]  <- daily_values_list[[i]]
#      daily_returns_mat[i, ] <- daily_returns_list[[i]]
#    }
#  
#    returns_xts_mat <- xts(t(daily_returns_mat), order.by = test_dates)
#    annual_returns <- as.numeric(Return.annualized(R = returns_xts_mat, scale = 252))
#    names(annual_returns) <- colnames(returns_xts_mat)
#    annual_vols <- as.numeric(StdDev.annualized(x = returns_xts_mat, scale = 252))
#    names(annual_vols) <- colnames(returns_xts_mat)
#    annual_sharp <- as.numeric(SharpeRatio.annualized(R = returns_xts_mat, scale = 252))
#    names(annual_sharp) <- colnames(returns_xts_mat)
#  
#    window_result_returns <- as.data.frame(t(annual_returns))
#    window_result_returns$Window <- window_index
#    rolling_annual_expected_returns <- rbind(rolling_annual_expected_returns, window_result_returns)
#  
#    window_result_sharpe <- as.data.frame(t(annual_sharp))
#    window_result_sharpe$Window <- window_index
#    rolling_annual_sharpe_ratios <- rbind(rolling_annual_sharpe_ratios, window_result_sharpe)
#  
#    window_result_vols <- as.data.frame(t(annual_vols))
#    window_result_vols$Window <- window_index
#    rolling_annual_volatilities <- rbind(rolling_annual_volatilities, window_result_vols)
#  
#    cat("Completed window", window_index,
#        ": Training rows [", train_start, "to", train_end,
#        "] (Dates:", format(train_data$Date[1], "%Y-%m-%d"),
#        "to", format(train_data$Date[nrow(train_data)], "%Y-%m-%d"),
#        "), Testing rows [", test_start, "to", test_end,
#        "] (Dates:", format(test_data$Date[1], "%Y-%m-%d"),
#        "to", format(test_data$Date[nrow(test_data)], "%Y-%m-%d"), ")\n")
#  }

