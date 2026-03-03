# Internal functions for Slab regression models: SR.
# These functions are not exported and are intended for internal use within the package's main function.

################################################################################
############################## Slab Regression #################################
################################################################################
# SR estimator
# Parameters:
#   v                - A scalar used to construct the vector u.
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda     - Regularized covariance matrix.
#   Sigma_lambda_inv - Inverse of Sigma_lambda.
# Returns:
#   A vector representing the SR estimator.
SR_ost <- function(v, est, sigma_square, Sigma_lambda, Sigma_lambda_inv) {
  u <- rep(v, length(est))
  a0_u <- sum(u^2)
  a1_u <- as.numeric(t(u) %*% (Sigma_lambda_inv %^% 1) %*% u)
  a2_u <- as.numeric(t(u) %*% (Sigma_lambda_inv %^% 2) %*% u)
  a3_u <- as.numeric(t(u) %*% (Sigma_lambda_inv %^% 3) %*% u)
  delta_u <- sigma_square * (a0_u * a3_u - a1_u * a2_u) + a3_u * (t(u) %*% est)^2
  mu_star <- (sigma_square * a2_u) / delta_u
  J <- matrix(1, length(est), length(est))
  if (delta_u > 0) {
    scalar_value <- as.numeric(mu_star / (1 + mu_star * a1_u))
    est_SR <- (diag(length(est)) - scalar_value * Sigma_lambda_inv %*% J) %*% est
  } else {
    est_SR <- (diag(length(est)) - Sigma_lambda_inv %*% J) %*% est
  }
  as.vector(est_SR)
}

# GSR estimator
# Parameters:
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda     - Regularized covariance matrix.
#   Sigma_lambda_inv - Inverse of Sigma_lambda (not used here but kept for consistency).
# Returns:
#   A vector representing the GSR estimator.
GSR_ost <- function(est, sigma_square, Sigma_lambda, Sigma_lambda_inv) {
  eigen_res <- eigen(Sigma_lambda, symmetric = TRUE)
  eigenvalues <- eigen_res$values
  U <- eigen_res$vectors
  beta_proj <- as.vector(crossprod(U, est))
  mu_l_star <- sigma_square / (beta_proj^2)
  scalar_coeff <- (mu_l_star / eigenvalues) / (1 + (mu_l_star / eigenvalues))
  adjustment_matrix <- diag(1, length(est)) - U %*% diag(scalar_coeff) %*% t(U)
  as.vector(adjustment_matrix %*% est)
}

