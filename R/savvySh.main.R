#' @title Slab and Shrinkage Linear Regression Estimation
#'
#' @description
#' This function estimates coefficients in a linear regression model using several shrinkage methods,
#' including Multiplicative Shrinkage, Slab Regression, Linear shrinkage, and Shrinkage Ridge Regression.
#' Each method gives estimators that balance bias and variance by applying shrinkage to the ordinary
#' least squares (OLS) solution. The shrinkage estimators are computed based on different assumptions
#' about the data.
#'
#' @usage savvySh(x, y, model_class = c("Multiplicative", "Slab", "Linear", "ShrinkageRR"),
#'                v = 1, lambda_vals = NULL, nlambda = 100, folds = 10,
#'                foldid = FALSE, include_Sh = FALSE, exclude = NULL)
#'
#' @param x A matrix of predictor variables.
#' @param y A vector of response variable.
#' @param model_class A character string specifying the shrinkage model to use. Options can choose from \code{"Multiplicative"},
#' \code{"Slab"}, \code{"Linear"}, and \code{"ShrinkageRR"}. The default is \code{"Multiplicative"}.
#'  If the user supplies more than one model, a warning is issued and only the first option is used.
#' @param v A numeric value controlling the strength of shrinkage for the \code{SR} estimator in the \code{"Slab"} model.
#' Must be a positive number. Default is 1.
#' @param lambda_vals A vector of \code{lambda} values for \code{RR}. This is used only when
#' multicollinearity (rank deficiency) is detected and \code{"ShrinkageRR"} is not selected.
#' If \code{NULL}, a default sequence is used.
#' @param nlambda The number of \code{lambda} values to use for cross-validation if \code{lambda_vals} is \code{NULL}.
#' Only used when multicollinearity is present and \code{"ShrinkageRR"} is not called. The default is \code{100}.
#' @param folds Number of folds for cross-validation in \code{RR}. This is applicable only if multicollinearity occurs
#' and \code{"ShrinkageRR"} is not chosen. The default is \code{10} and must be an integer \code{>= 3}.
#' @param foldid Logical. If \code{TRUE}, saves the fold assignments in the output when multicollinearity is detected
#' and \code{"ShrinkageRR"} is not used. The default is \code{FALSE}.
#' @param include_Sh Logical. If \code{TRUE}, includes the Sh estimator in the \code{"Multiplicative"} model. The default is \code{FALSE}.
#' @param exclude A vector specifying columns to exclude from the predictors. The default is \code{NULL}.
#'
#' @details
#' The \emph{Slab and Shrinkage Linear Regression Estimation} methodology provides four classes of shrinkage estimators
#' that reduce variance in the OLS solution by introducing a small, structured bias. These methods handle overfitting,
#' collinearity, and high-dimensional scenarios by controlling how and where the coefficients are shrunk. Each class offers a distinct strategy
#' for controlling instability and improving mean squared error (MSE) in linear models, tailored for different modeling contexts specified
#' in the \code{model_class} argument. Note that if the user provides more than one option in \code{model_class}, only the first option is used,
#' and a warning is issued.
#'
#' \strong{Model Classes:}
#' \describe{
#'
#'   \item{\strong{Multiplicative Shrinkage:}}{
#'     This class includes three estimators that use the \code{OLS} coefficients as a starting point and apply
#'     \emph{multiplicative} adjustments:
#'
#'     \describe{
#'       \item{\code{St - }}{\emph{Stein estimator}, which shrinks all coefficients toward zero by a single global factor.
#'         This aims to reduce MSE while keeping the overall bias fairly uniform across coefficients.}
#'       \item{\code{DSh - }}{\emph{Diagonal Shrinkage}, assigning an individual factor to each coefficient based on its
#'         variance. This yields more targeted shrinkage than the global approach and often achieves a lower MSE.}
#'       \item{\code{Sh - }}{\emph{Shrinkage estimator} that solves a Sylvester equation for a full (non-diagonal) shrinkage matrix.
#'         It is more flexible but also more computationally demanding. Included only if \code{include_Sh = TRUE}.}
#'     }
#'   }
#'
#'   \item{\strong{Slab Regression:}}{
#'     \emph{Slab Regression} applies an adaptive quadratic penalty term to the \code{OLS} objective:
#'
#'     \describe{
#'       \item{\code{SR - }}{\emph{Simple Slab Regression}, which modifies the \code{OLS} objective by
#'         adding a penalty in a fixed direction (often the constant vector). This penalty is controlled by \code{v}
#'         and does not require cross-validation. It can be viewed as a special case of the \code{generalized lasso}
#'         but focuses on smooth (quadratic) rather than \eqn{\ell_1}{L1} regularization.}
#'       \item{\code{GSR - }}{\emph{Generalized Slab Regression}, extending SR by allowing shrinkage along multiple directions.
#'         Typically, these directions correspond to the eigenvectors of the design covariance matrix, effectively shrinking
#'         principal components.}
#'     }
#'   }
#'
#'   \item{\strong{Linear Shrinkage:}}{
#'     The \emph{Linear Shrinkage (LSh)} estimator forms a convex combination of the OLS estimator (through the origin)
#'     and a target estimator that assumes uncorrelated predictors (\emph{diagonal} approximation of the covariance).
#'     This approach is simpler than a full matrix method and is well-suited for standardized data where the intercept
#'     is not needed.
#'   }
#'
#'   \item{\strong{Shrinkage Ridge Regression:}}{
#'     The \emph{Shrinkage Ridge Regression (SRR)} extends standard \code{RR} by shrinking the design covariance
#'     matrix toward a spherical target (i.e., a diagonal matrix with equal entries). This additional regularization
#'     stabilizes the eigenvalues and yields more robust coefficient estimates, particularly when the predictors lie
#'     close to a low-dimensional subspace.
#'   }
#' }
#'
#' @return A list containing the following elements:
#' \item{call}{The matched function call.}
#' \item{model}{The data frame of \code{y} and \code{x} used in the analysis.}
#' \item{\code{optimal_lambda}}{If \code{x} is full rank, this value is \code{0}.
#' If \code{x} is rank-deficient, it is the chosen \code{RR} \code{lambda} from cross-validation.}
#' \item{model_class}{The selected model class.}
#' \item{coefficients}{A list of estimated coefficients for each applicable estimator in the \code{model_class}.}
#' \item{fitted_values}{A list of fitted values for each estimator.}
#' \item{pred_MSE}{A list of prediction MSEs for each estimator.}
#' \item{ridge_results (optional)}{
#'  A list containing detailed results from \code{RR}, used when multicollinearity (rank deficiency)
#'  is detected in \code{x} and the \code{"ShrinkageRR"} is not called. This element is included only when \code{RR} is applied instead of \code{OLS}
#'  due to the rank deficiency of \code{x}. It contains:
#'  \describe{
#'  \item{\code{lambda_range}}{The range of \code{lambda} values used in the \code{RR} cross-validation.}
#'  \item{\code{cvm}}{A vector of cross-validated MSEs for each \code{lambda} in \code{lambda_range}.}
#'  \item{\code{cvsd}}{The standard deviation of the cross-validated MSEs for each \code{lambda}.}
#'  \item{\code{ridge_coefficients}}{A matrix of coefficients from \code{RR} at each \code{lambda} value,
#'  with each column representing the coefficients corresponding to a specific \code{lambda}.}
#'  }
#' }
#'
#' @examples
#' # 1. Simple Multiplicative Shrinkage example
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- rnorm(100)
#' fit_mult <- savvySh(x, y, model_class = "Multiplicative", include_Sh = TRUE)
#' print(fit_mult)
#'
#' # 2. Slab Regression example
#' fit_slab <- savvySh(x, y, model_class = "Slab", v = 2)
#' coef(fit_slab, estimator = "GSR")
#'
#' # 3. Linear Shrinkage (standardized data recommended)
#' x_centered <- scale(x, center = TRUE, scale = FALSE)
#' y_centered <- scale(y, center = TRUE, scale = FALSE)
#' fit_linear <- savvySh(x_centered, y_centered, model_class = "Linear")
#'
#' # 4. Shrinkage Ridge Regression
#' fit_srr <- savvySh(x, y, model_class = "ShrinkageRR")
#' predict(fit_srr, newx = matrix(rnorm(10 * 5), 10, 5), type = "response")
#'
#' @author
#' Ziwei Chen, Vali Asimit, Marina Anca Cidota, Jennifer Asimit\cr
#' Maintainer: Ziwei Chen <ziwei.chen.3@citystgeorges.ac.uk>
#'
#' @references
#' Asimit, V., Cidota, M. A., Chen, Z., & Asimit, J. (2025).
#' \emph{Slab and Shrinkage Linear Regression Estimation}.
#' Retrieved from \url{https://openaccess.city.ac.uk/id/eprint/35005/}
#'
#' @importFrom Matrix rankMatrix Schur
#' @importFrom mnormt pd.solve
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom MASS mvrnorm ginv
#' @importFrom expm "%^%"
#' @importFrom stats coef lm predict optimize
#' @export
savvySh <- function(x, y, model_class = c("Multiplicative", "Slab", "Linear", "ShrinkageRR"),
                    v = 1, lambda_vals = NULL, nlambda = 100, folds = 10, foldid = FALSE,
                    include_Sh = FALSE, exclude = NULL) {

  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) {
    stop("The number of rows in x must match the length of y.")
  }
  if (anyNA(x) || anyNA(y)) {
    stop("x or y contains missing values. Please impute or handle missing data before proceeding.")
  }
  valid_models <- c("Multiplicative", "Slab", "Linear", "ShrinkageRR")
  if (length(model_class) > 1) {
    warning("model_class should be a single option. Using the first element: ", model_class[1])
    model_class <- model_class[1]
  } else {
    model_class <- match.arg(model_class, choices = model_class)
  }
  if (!model_class %in% valid_models) {
    stop(paste("Invalid model_class specified. Choose from:", paste(valid_models, collapse = ", ")))
  }
  if (include_Sh && model_class != "Multiplicative") {
    warning("include_Sh is only applicable when model_class is 'Multiplicative'. Setting include_Sh to FALSE.")
    include_Sh <- FALSE
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste("V", seq_len(ncol(x)), sep = "")
  }
  if (!is.null(exclude)) {
    if (anyNA(exclude) || any(exclude == 0)) {
      stop("Exclusion must not contain NA or zero values.")
    } else if (!is.vector(exclude) || is.matrix(exclude) || is.list(exclude)) {
      stop("Exclusion must be a vector.")
    } else if (is.numeric(exclude)) {
      if (any(exclude > ncol(x))) {
        stop("Exclusion indices are out of bounds.")
      }
      x <- x[, -exclude, drop = FALSE]
    } else if (is.character(exclude)) {
      if (!all(exclude %in% colnames(x))) {
        stop("Exclusion names must match column names.")
      }
      x <- x[, !colnames(x) %in% exclude, drop = FALSE]
    }
  }

  nobs <- nrow(x)
  nvars <- ncol(x)
  x_tilde <- cbind(1, x)
  if (nvars >= nobs) {
    warning("Number of features in x is less than the number of observations.")
  }
  full_rank <- Matrix::rankMatrix(x)[1] == nvars
  if (!full_rank && model_class != "ShrinkageRR") {
    if (is.null(lambda_vals)) {
      lambda_vals <- 10^seq(-6, 2, length.out = nlambda)
    } else if (length(lambda_vals) < 2) {
      stop("Need more than one value of lambda_vals for cross-validation.")
    }
    if (!is.numeric(folds) || folds < 3 || round(folds) != folds) {
      stop("Number of folds must be an integer greater than or equal to 3.")
    }
    fold_ids <- sample(rep(1:folds, length.out = nobs))
  }

  if (model_class %in% c("Multiplicative", "Slab")) {
    if (!full_rank) {
      warning("Multicollinearity detected: switched to RR instead of unbiased OLS estimation.")
      glmnet_fit <- cv.glmnet(x, y, alpha = 0, lambda = lambda_vals, nfolds = folds)
      optimal_lambda <- glmnet_fit$lambda.min
      est <- as.vector(coef(glmnet_fit, s = "lambda.min"))
      fitted_values <- as.vector(predict(glmnet_fit, newx = x, s = optimal_lambda))
      df_eff <- sum((svd(x)$d)^2 / ((svd(x)$d)^2 + optimal_lambda))
      RSS <- sum((y - fitted_values)^2)
      sigma_square <- (RSS / (nobs - df_eff))^2
      ridge_results <- list(
        lambda_range = glmnet_fit$lambda,
        cvm = glmnet_fit$cvm,
        cvsd = glmnet_fit$cvsd,
        ridge_coefficients = as.matrix(coef(glmnet_fit, s = glmnet_fit$lambda))
      )
    } else {
      lm_fit <- lm(y ~ x)
      est <- as.vector(lm_fit$coefficients)
      sigma_square <- (summary(lm_fit)$sigma)^2
      optimal_lambda <- 0
    }
    Sigma_lambda <- t(x_tilde) %*% x_tilde + optimal_lambda * diag(ncol(x_tilde))
    Sigma_lambda <- (Sigma_lambda + t(Sigma_lambda)) / 2
    Sigma_lambda_inv <- pd.solve(Sigma_lambda)

    if (model_class == "Multiplicative") {
      est_St <- St_ost(est, sigma_square, Sigma_lambda_inv)
      St_fitted_values <- as.vector(x_tilde %*% est_St)
      pred_MSE_St <- mean((y - St_fitted_values)^2)

      est_DSh <- DSh_ost(est, sigma_square, Sigma_lambda_inv)
      DSh_fitted_values <- as.vector(x_tilde %*% est_DSh)
      pred_MSE_DSh <- mean((y - DSh_fitted_values)^2)

      if (include_Sh) {
        est_Sh <- Sh_ost(est, sigma_square, Sigma_lambda_inv)
        Sh_fitted_values <- as.vector(x_tilde %*% est_Sh)
        pred_MSE_Sh <- mean((y - Sh_fitted_values)^2)
      }

    } else if (model_class == "Slab") {
      if (!is.numeric(v) || v <= 0) {
        stop("Error: 'v' must be a positive number. This parameter is required only when using the SR estimator.")
      }
      est_SR <- SR_ost(v, est, sigma_square, Sigma_lambda, Sigma_lambda_inv)
      SR_fitted_values <- as.vector(x_tilde %*% est_SR)
      pred_MSE_SR <- mean((y - SR_fitted_values)^2)

      est_GSR <- GSR_ost(est, sigma_square, Sigma_lambda, Sigma_lambda_inv)
      GSR_fitted_values <- as.vector(x_tilde %*% est_GSR)
      pred_MSE_GSR <- mean((y - GSR_fitted_values)^2)
    }

  } else if (model_class == "Linear") {
    if (nvars == 1) {
      stop("The number of covariates must be at least two for this model.")
    }
    centered_x <- scale(x, center = TRUE, scale = FALSE)
    centered_y <- scale(y, center = TRUE, scale = FALSE)

    if (!full_rank) {
      warning("Multicollinearity detected: switched to RR on centered data instead of unbiased OLS estimation.")
      glmnet_fit <- cv.glmnet(centered_x, centered_y, alpha = 0, lambda = lambda_vals, nfolds = folds)
      optimal_lambda <- glmnet_fit$lambda.min
      est <- as.vector(coef(glmnet_fit, s = "lambda.min"))[-1]
      fitted_values <- as.vector(predict(glmnet_fit, newx = centered_x, s = optimal_lambda))
      df_eff <- sum((svd(centered_x)$d)^2 / ((svd(centered_x)$d)^2 + optimal_lambda))
      RSS <- sum((centered_y - fitted_values)^2)
      sigma_square <- (RSS / (nobs - df_eff))^2
      ridge_results <- list(
        lambda_range = glmnet_fit$lambda,
        cvm = glmnet_fit$cvm,
        cvsd = glmnet_fit$cvsd,
        ridge_coefficients = as.matrix(coef(glmnet_fit, s = glmnet_fit$lambda))
      )
    } else {
      lm_fit <- lm(centered_y ~ centered_x - 1)
      est <- as.vector(lm_fit$coefficients)
      sigma_square <- (summary(lm_fit)$sigma)^2
      optimal_lambda <- 0
    }
    Sigma_lambda <- t(centered_x) %*% centered_x + optimal_lambda * diag(ncol(centered_x))
    Sigma_lambda <- (Sigma_lambda + t(Sigma_lambda)) / 2
    Sigma_lambda_inv <- pd.solve(Sigma_lambda)
    Sigma_tilde <- diag(diag(Sigma_lambda))

    est_LSh <- LSh_ost(est, sigma_square, Sigma_lambda, Sigma_lambda_inv, Sigma_tilde)
    LSh_fitted_values <- as.vector(centered_x %*% est_LSh)
    pred_MSE_LSh <- mean((centered_y - LSh_fitted_values)^2)

  } else if (model_class == "ShrinkageRR") {
    Sigma_lambda <- t(x_tilde) %*% x_tilde
    est_SRR <- SRR_ost(x_tilde, y, Sigma_lambda)
    SRR_fitted_values <- as.vector(x_tilde %*% est_SRR)
    pred_MSE_SRR <- mean((y - SRR_fitted_values)^2)
  }

  results <- list(
    call = match.call(),
    model = data.frame(y, x),
    model_class = model_class
  )

  if (model_class == "Multiplicative") {
    results$optimal_lambda <- optimal_lambda
    results$coefficients <- list(St = est_St, DSh = est_DSh)
    results$fitted_values <- list(St = St_fitted_values, DSh = DSh_fitted_values)
    results$pred_MSE <- list(St = pred_MSE_St, DSh = pred_MSE_DSh)
    if (include_Sh) {
      results$coefficients$Sh <- est_Sh
      results$fitted_values$Sh <- Sh_fitted_values
      results$pred_MSE$Sh <- pred_MSE_Sh
    }

  } else if (model_class == "Slab") {
    results$optimal_lambda <- optimal_lambda
    results$coefficients <- list(SR = est_SR, GSR = est_GSR)
    results$fitted_values <- list(SR = SR_fitted_values, GSR = GSR_fitted_values)
    results$pred_MSE <- list(SR = pred_MSE_SR, GSR = pred_MSE_GSR)

  } else if (model_class == "Linear") {
    results$optimal_lambda <- optimal_lambda
    results$coefficients <- list(LSh = est_LSh)
    results$fitted_values <- list(LSh = LSh_fitted_values)
    results$pred_MSE <- list(LSh = pred_MSE_LSh)

  } else if (model_class == "ShrinkageRR") {
    results$coefficients <- list(SRR = est_SRR)
    results$fitted_values <- list(SRR = SRR_fitted_values)
    results$pred_MSE <- list(SRR = pred_MSE_SRR)
  }

  if (!full_rank && model_class != "ShrinkageRR") {
    results$ridge_results <- ridge_results
    if (foldid) {
      results$fold_assignments <- fold_ids
    }
  }

  class(results) <- "savvySh_model"
  return(results)
}


