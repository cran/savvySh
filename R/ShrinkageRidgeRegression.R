# Internal functions for Shrinkage ridge shrinkage models: SRR.
# These functions are not exported and are intended for internal use within the package's main function.

################################################################################
########################## Shrinkage Ridge Regression ##########################
################################################################################
# Helper function to compute H(rho) for the SRR model.
# Parameters:
#   rho       - Shrinkage parameter (0 <= rho <= 1).
#   nobs      - Number of observations.
#   p_plus_1  - Number of predictors plus one for the intercept.
#   yty       - Squared norm of the response vector (y^T * y).
#   lambda    - Vector of eigenvalues of the covariance matrix.
#   V_k       - Variance terms associated with each eigenvalue.
#   V_star    - Average eigenvalue of the covariance matrix.
# Returns:
#   A numeric value representing H(rho).
H_function <- function(rho, nobs, p_plus_1, yty, lambda, V_k, V_star) {
  denom  <- (1 - rho) * lambda + rho * V_star
  denom2 <- denom^2
  denom3 <- denom^3
  term1 <- sum(lambda / denom2)
  term2 <- sum(V_k / denom)
  term3 <- sum(lambda * V_k / denom2)
  term4 <- sum(V_k * rho^2 * (lambda - V_star)^2 / denom3)
  H_rho <- (1 / (nobs - p_plus_1)) * (yty - 2 * term2 + term3) * term1 + term4
  H_rho
}

# Function to find the optimal rho for the SRR model.
# Parameters:
#   nobs      - Number of observations.
#   p_plus_1  - Number of predictors plus one for the intercept.
#   yty       - Squared norm of the response vector (y^T * y).
#   lambda    - Vector of eigenvalues of the covariance matrix.
#   V_k       - Variance terms associated with each eigenvalue.
#   V_star    - Average eigenvalue of the covariance matrix.
# Returns:
#   A numeric value representing the optimal rho that minimizes H(rho).
find_rho_star <- function(nobs, p_plus_1, yty, lambda, V_k, V_star) {
  c_function <- function(rho) H_function(rho, nobs, p_plus_1, yty, lambda, V_k, V_star)
  rho_opt <- optimize(c_function, interval = c(0, 1))
  rho_opt$minimum
}

# Main SRR function.
# Computes the Shrinkage Ridge Regression (SRR) estimator.
# Parameters:
#   x_tilde   - Design matrix including the intercept column.
#   y         - Response vector.
#   Sigma_lambda - Regularized covariance matrix (X^T X + lambda * I).
# Returns:
#   A vector representing the SRR estimator (estimated coefficients).
SRR_ost <- function(x_tilde, y, Sigma_lambda) {
  p_plus_1 <- ncol(x_tilde)
  eigen_res <- eigen(Sigma_lambda, symmetric = TRUE)
  lambda <- eigen_res$values
  mu     <- eigen_res$vectors
  V_star <- (1 / p_plus_1) * sum(lambda)
  yty <- as.numeric(crossprod(y))
  proj <- as.vector(t(y) %*% (x_tilde %*% mu))
  V_k <- proj^2
  rho_star <- find_rho_star(nrow(x_tilde), p_plus_1, yty, lambda, V_k, V_star)
  diag_entries <- 1 / ((1 - rho_star) * lambda + rho_star * V_star)
  beta_SRR <- mu %*% diag(diag_entries, nrow = length(diag_entries)) %*% t(mu)
  est_SRR <- beta_SRR %*% t(x_tilde) %*% y
  as.vector(est_SRR)
}


