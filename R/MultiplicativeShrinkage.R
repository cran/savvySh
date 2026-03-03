# Internal functions for multiplicative shrinkage models: St, DSh, and Sh.
# These functions are not exported and are intended for internal use within the package's main function.

################################################################################
########################## Multiplicative Shrinkage ############################
################################################################################
# Main St function.
# Computes the Stein estimator (St) using multiplicative shrinkage.
# Parameters:
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda_inv - Inverse of the regularized covariance matrix.
# Returns:
#   A vector representing the St estimator.
St_ost <- function(est, sigma_square, Sigma_lambda_inv) {
  M0 <- sigma_square * sum(diag(Sigma_lambda_inv))
  a_star <- sum(est^2) / (sum(est^2) + M0)
  as.vector(a_star * est)
}

# Main DSh function.
# Computes the Diagonal Shrinkage (DSh) estimator using multiplicative shrinkage.
# Parameters:
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda_inv - Inverse of the regularized covariance matrix.
# Returns:
#   A vector representing the DSh estimator.
DSh_ost <- function(est, sigma_square, Sigma_lambda_inv) {
  b_star <- (est^2) / (est^2 + sigma_square * diag(Sigma_lambda_inv))
  as.vector(b_star * est)
}

# Sylvester equation solver.
# Solves the Sylvester equation A * X + X * B = C using a Schur decomposition method.
# This code for solving the Sylvester equation was adapted from the biADMM R code,
# available at: https://github.com/sakuramomo1005/biADMM/blob/master/R/sylvester.R
# Original Author: sakuramomo1005
# Parameters:
#   A      - Numeric matrix.
#   B      - Numeric matrix.
#   C      - Numeric matrix.
#   tol    - Tolerance level for small off-diagonal entries (default = 1e-4).
# Returns:
#   A matrix representing the solution X.
sylvester <- function(A, B, C, tol = 1e-4) {
  A1 <- Schur(A)
  Q1 <- A1$Q
  R1 <- A1$T

  A2 <- Schur(B)
  Q2 <- A2$Q
  R2 <- A2$T

  C_trans <- t(Q1) %*% C %*% Q2
  n <- nrow(R2)
  p <- ncol(R1)
  X <- matrix(0, p, n)
  I <- diag(p)
  k <- 1

  while (k <= n) {
    if (k < n && abs(R2[k+1, k]) >= tol) {
      # Process 2x2 block.
      r11 <- R2[k, k]
      r12 <- R2[k, k+1]
      r21 <- R2[k+1, k]
      r22 <- R2[k+1, k+1]
      if (k == 1) {
        temp <- matrix(0, p, 2)
      } else {
        temp <- X[, 1:(k-1), drop = FALSE] %*% R2[1:(k-1), k:(k+1), drop = FALSE]
      }
      b_block <- C_trans[, k:(k+1), drop = FALSE] - temp
      A_mat <- R1 %*% R1 + (r11 + r22) * R1 + (r11 * r22 - r12 * r21) * I
      X[, k:(k+1)] <- ginv(A_mat) %*% b_block
      k <- k + 2
    } else {
      # Process single column.
      if (k == 1) {
        temp <- matrix(0, p, 1)
      } else {
        temp <- X[, 1:(k-1), drop = FALSE] %*% R2[1:(k-1), k, drop = FALSE]
      }
      b_col <- C_trans[, k, drop = FALSE] - temp
      X[, k] <- ginv(R1 + R2[k, k] * I) %*% b_col
      k <- k + 1
    }
  }
  Q1 %*% X %*% t(Q2)
}

# Main Sh function.
# Computes the Shrinkage (Sh) estimator by solving a Sylvester equation to obtain
# a targeted shrinkage solution.
# Parameters:
#   est              - Initial coefficient estimate vector.
#   sigma_square     - Variance of noise.
#   Sigma_lambda_inv - Inverse of the regularized covariance matrix.
# Returns:
#   A vector representing the Sh estimator.
Sh_ost <- function(est, sigma_square, Sigma_lambda_inv) {
  B <- est %*% t(est)
  C_star <- sylvester(Sigma_lambda_inv, B, B)
  as.vector(C_star %*% est)
}


