# Internal functions for Linear shrinkage models: LSh.
# These functions are not exported and are intended for internal use within the package's main function.

################################################################################
############################# Linear Shrinkage  ################################
################################################################################
# LSh estimator
# Parameters:
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda     - Regularized covariance matrix.
#   Sigma_lambda_inv - Inverse of Sigma_lambda.
#   Sigma_tilde      - Shrinkage target covariance matrix.
# Returns:
#   A vector representing the LSh estimator.
LSh_ost <- function(est, sigma_square, Sigma_lambda, Sigma_lambda_inv, Sigma_tilde) {
  I_p <- diag(length(est))
  Sigma_tilde_inv <- pd.solve(Sigma_tilde)
  t1 <- sigma_square * sum(diag(Sigma_tilde_inv))
  t2 <- sigma_square * sum(diag(Sigma_lambda_inv))
  diff_matrix <- Sigma_tilde_inv %*% Sigma_lambda - I_p
  t3 <- as.numeric(t(est) %*% (diff_matrix %*% diff_matrix) %*% est)
  rho_star <- (t2 - t1) / (t2 - t1 + t3)
  Sigma_rho_star <- rho_star * Sigma_tilde_inv %*% Sigma_lambda + (1 - rho_star) * I_p
  as.vector(Sigma_rho_star %*% est)
}

