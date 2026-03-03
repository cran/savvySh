#' @title Print a Slab and Shrinkage Model Summary
#'
#' @description
#' Displays a concise summary of a fitted \code{savvySh_model} object, including the original
#' function call, the chosen model class, the number of non-zero coefficients per estimator,
#' and the optimal \code{lambda} value (if applicable). Additionally, it prints the coefficients for
#' the specified estimator(s) with user-specified precision.
#'
#' @param x A fitted \code{savvySh_model} object returned by \code{savvySh}.
#' @param digits An integer specifying the number of significant digits to display when printing
#'   coefficients and \code{lambda}. Defaults to \code{max(3, getOption("digits") - 3)}.
#' @param estimator A character vector naming one or more estimators for which coefficients should be printed.
#'   Valid names are those present in \code{x$coefficients} (e.g., \code{"St"}, \code{"DSh"}, \code{"Sh"},
#'   \code{"SR"}, \code{"GSR"}, or \code{"ShrinkageRR"}). If \code{NULL}, coefficients for all estimators are printed.
#' @param ... Additional arguments passed to \code{\link[base]{print}} (currently unused).
#'
#' @details
#' This print method provides a quick diagnostic of the fitted model by showing:
#' \describe{
#'   \item{Summary Metrics}{A table that includes, for each estimator, the number of non-zero coefficients
#'       and the optimal \code{lambda} (if applicable).}
#'   \item{Coefficients}{For each selected estimator, the coefficients are printed with appropriate names:
#'       if an intercept is present, it is labeled \code{(Intercept)} and the remaining
#'       coefficients are labeled according to the predictor names.}
#' }
#' If the user does not specify an estimator using the \code{estimator} argument, the function prints
#' information for all available estimators stored in the model. If one or more estimators are specified,
#' only those are printed, after verifying that they exist in \code{x$coefficients}.
#'
#' The method invisibly returns a summary \code{data.frame} containing key metrics for each estimator.
#'
#' @return
#' Invisibly returns a \code{data.frame} summarizing each selected estimator's name, number of non-zero
#' coefficients, and the final \code{optimal_lambda} (if any).
#'
#' @seealso
#'   \code{\link{savvySh}} for fitting slab and shrinkage linear models,
#'   \code{\link{coef.savvySh_model}} and \code{\link{predict.savvySh_model}} for extracting coefficients
#'   and generating predictions.
#'
#' @examples
#' # Generate simulated data for demonstration
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- rnorm(100)
#'
#' # Fit a Multiplicative shrinkage model
#' fit <- savvySh(x, y, model_class = "Multiplicative", include_Sh = TRUE)
#'
#' # Default print: shows summary metrics and coefficients for all estimators
#' print(fit)
#'
#' # Print with specific digits and only for one estimator
#' print(fit, digits = 4, estimator = "St")
#'
#' @author
#' Ziwei Chen, Vali Asimit, Marina Anca Cidota, Jennifer Asimit\cr
#' Maintainer: Ziwei Chen <ziwei.chen.3@citystgeorges.ac.uk>
#'
#' @method print savvySh_model
#' @export
print.savvySh_model <- function(x, digits = max(3, getOption("digits") - 3), estimator = NULL, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Model Class: ", x$model_class, "\n\n")

  if (is.null(x$coefficients) || !is.list(x$coefficients)) {
    stop("The 'savvySh_model' object does not contain valid coefficients.")
  }
  if (is.null(estimator)) {
    estimator <- names(x$coefficients)
  } else {
    if (!is.character(estimator)) {
      stop("Estimator must be a character vector.")
    }
    if (!all(estimator %in% names(x$coefficients))) {
      stop("The specified estimator(s) '", paste(estimator, collapse = ", "),
           "' are not available. Choose from: ", paste(names(x$coefficients), collapse = ", "))
    }
  }

  summary_list <- lapply(estimator, function(est) {
    coefs <- x$coefficients[[est]]
    list(
      "Estimator" = est,
      "Non-Zero Coefficients" = sum(coefs != 0),
      "Optimal Lambda" = ifelse(is.null(x$optimal_lambda), NA, signif(x$optimal_lambda, digits))
    )
  })
  summary_output <- do.call(rbind, lapply(summary_list, as.data.frame))
  print(summary_output, row.names = FALSE)

  cat("\nCoefficients:\n")
  first_coefs <- x$coefficients[[estimator[1]]]
  intercept_present <- (ncol(x$model) == length(first_coefs))

  for (est in estimator) {
    cat("\nEstimator: ", est, "\n")
    coefs <- x$coefficients[[est]]
    coef_names <- if (intercept_present) {
      c("(Intercept)", colnames(x$model)[-1])
    } else {
      paste0("V", seq_along(coefs))
    }
    names(coefs) <- coef_names
    coef_df <- data.frame(
      Coefficient = coef_names,
      Estimate = signif(coefs, digits)
    )
    print(coef_df, row.names = FALSE)
  }

  invisible(summary_output)
}
