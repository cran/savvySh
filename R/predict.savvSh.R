#' @title Predict Method for Slab and Shrinkage Linear Regression Models
#'
#' @description
#' Generate predictions (fitted values) or extract regression coefficients from a
#' \code{savvySh_model} object returned by \code{savvySh}. This function allows you to
#' specify one or more shrinkage estimators (via the \code{estimator} parameter) available
#' in the model. If no estimator is specified, all available estimators are used and their
#' results are returned in a named list.
#'
#' @param object A fitted \code{savvySh_model} object produced by \code{savvySh}.
#' @param newx A numeric matrix of new predictor data for which to generate predictions.
#'   This argument is required if \code{type = "response"} and is ignored if \code{type = "coefficients"}.
#' @param type A character string specifying the output type. Options are \code{"response"} to return
#'   predicted values and \code{"coefficients"} to extract regression coefficient vectors. Defaults to \code{"response"}.
#' @param estimator A character vector naming one or more shrinkage estimator(s) to use.
#'   These must match names present in \code{object$coefficients}. If \code{NULL},
#'   all available estimators are used.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The behavior depends on the value of \code{type}:
#' \describe{
#'   \item{\code{"response"}:}{Generates predicted values using the coefficient estimates
#'       from the specified shrinkage estimator(s) for new data supplied via \code{newx}.}
#'   \item{\code{"coefficients"}:}{Extracts the regression coefficient vector(s) corresponding
#'       to the specified estimator(s). Coefficient names are assigned based on whether an intercept
#'       is present (for Linear shrinkage, no intercept).}
#' }
#'
#' If no \code{estimator} is specified, the function returns results for all available estimators
#' as a named list. If a single estimator is specified (or only one is provided in the vector), the result
#' is returned as a numeric vector (for coefficients) or a numeric vector of predictions (for response).
#'
#' @return
#' If \code{type = "response"}, the function returns:
#' \itemize{
#'   \item A numeric vector of predicted values if exactly one estimator is specified;
#'   \item Otherwise, a named list of numeric vectors, one for each specified estimator.
#' }
#' If \code{type = "coefficients"}, the function returns:
#' \itemize{
#'   \item A named numeric vector of regression coefficients if exactly one estimator is specified;
#'   \item Otherwise, a named list of numeric vectors corresponding to each specified estimator.
#' }
#'
#' @seealso
#'   \code{\link{savvySh}} for fitting slab and shrinkage linear models,
#'   \code{\link{coef.savvySh_model}} for direct coefficient extraction.
#'
#' @examples
#' # Generate simulated data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- rnorm(100)
#'
#' # Fit a Multiplicative shrinkage model
#' fit <- savvySh(x, y, model_class = "Multiplicative")
#'
#' # Generate predictions for new data
#' new_x <- matrix(rnorm(10 * 5), 10, 5)
#' preds <- predict(fit, newx = new_x, type = "response")
#'
#' # Extract coefficients for specific estimators
#' coefs_st <- predict(fit, type = "coefficients", estimator = "St")
#'
#' @author
#' Ziwei Chen, Vali Asimit, Marina Anca Cidota, Jennifer Asimit\cr
#' Maintainer: Ziwei Chen <ziwei.chen.3@citystgeorges.ac.uk>
#'
#' @method predict savvySh_model
#' @export
predict.savvySh_model <- function(object, newx = NULL, type = c("response", "coefficients"),
                                  estimator = NULL, ...) {
  valid_types <- c("response", "coefficients")
  if (length(type) != 1) {
    stop("Please specify exactly one type: 'response' or 'coefficients'.")
  }
  if (!type %in% valid_types) {
    stop("Invalid type specified. Use 'response' or 'coefficients'.")
  }
  if (is.null(object$coefficients) || !is.list(object$coefficients)) {
    stop("Invalid 'savvySh_model' object: coefficients are missing or improperly structured.")
  }

  if (is.null(estimator)) {
    if (length(names(object$coefficients)) > 1) {
      estimator <- names(object$coefficients)
      warning("No estimator specified. Returning results for all available estimators: ",
              paste(estimator, collapse = ", "))
    } else {
      estimator <- names(object$coefficients)
    }
  } else {
    if (!is.character(estimator)) {
      stop("Estimator must be a character vector of estimator names.")
    }
    if (!all(estimator %in% names(object$coefficients))) {
      stop("The specified estimator(s) '", paste(estimator, collapse = ", "),
           "' are not available. Choose from: ", paste(names(object$coefficients), collapse = ", "))
    }
  }

  first_coefs <- object$coefficients[[estimator[1]]]
  intercept_present <- (ncol(object$model) == length(first_coefs))

  if (type == "coefficients") {
    coefs_list <- lapply(estimator, function(est_name) {
      coefs <- object$coefficients[[est_name]]
      if (is.null(coefs) || length(coefs) == 0) {
        stop(paste("The coefficients for the estimator", est_name, "are missing or empty."))
      }
      coef_names <- if (intercept_present) {
        c("(Intercept)", paste0("V", seq_len(ncol(object$model) - 1)))
      } else {
        paste0("V", seq_len(length(coefs)))
      }
      names(coefs) <- coef_names
      coefs
    })
    if (length(estimator) == 1) {
      return(coefs_list[[1]])
    } else {
      names(coefs_list) <- estimator
      return(coefs_list)
    }
  }

  if (type == "response") {
    if (is.null(newx)) {
      stop("You need to supply a matrix of new values for 'newx'.")
    }
    newx <- as.matrix(newx)
    nvars <- if (intercept_present) length(first_coefs) - 1 else length(first_coefs)
    if (ncol(newx) != nvars) {
      stop("The number of columns in 'newx' must match the number of predictors in the model.")
    }
    predictions <- lapply(estimator, function(est_name) {
      coefs <- object$coefficients[[est_name]]
      inter <- if (intercept_present) coefs[1] else 0
      beta <- if (intercept_present) coefs[-1] else coefs
      as.vector(inter + newx %*% beta)
    })
    if (length(estimator) == 1) {
      return(predictions[[1]])
    } else {
      names(predictions) <- estimator
      return(predictions)
    }
  }
}
