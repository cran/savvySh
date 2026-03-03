#' @title Summarize a Slab and Shrinkage Linear Regression Model
#'
#' @description
#' Provides a comprehensive summary for one or more shrinkage estimators contained within a
#' \code{savvySh_model} object produced by \code{savvySh}. The summary includes estimated coefficients,
#' confidence intervals, residual statistics, R-squared measures, F-statistics, and information criteria (AIC, BIC)
#' for each specified estimator.
#'
#' @param object A fitted model object of class \code{savvySh_model}, produced by \code{savvySh}.
#' @param estimator A character vector naming one or more estimators to summarize (e.g., \code{"St"}, \code{"DSh"},
#'   \code{"SR"}, \code{"GSR"}, \code{"Sh"}, etc.). If \code{NULL} (default), summaries for all available estimators are printed.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' For each estimator present in the \code{savvySh_model} object (or for the user-specified subset), this function computes:
#' \itemize{
#'   \item A summary of the residual distribution (quantiles).
#'   \item A coefficient table including estimates, standard errors, t-values, p-values, and confidence intervals.
#'   \item Residual standard error and degrees of freedom.
#'   \item R-squared and adjusted R-squared measures.
#'   \item F-statistic (and its p-value) for testing overall regression significance.
#'   \item Information criteria (AIC, BIC) and deviance for model fit.
#' }
#' These results are printed in sequence for the selected estimator(s). If no estimator is specified,
#' summaries for all available estimators are printed.
#'
#' @return
#' Invisibly returns a \code{data.frame} summarizing key metrics for each estimator (including estimator name,
#' number of non-zero coefficients, and optimal \code{lambda} if available).
#'
#' @seealso
#'   \code{\link{savvySh}} for fitting slab and shrinkage linear models,
#'   \code{\link{predict.savvySh_model}} for generating predictions,
#'   \code{\link{coef.savvySh_model}} for extracting coefficients directly.
#'
#' @examples
#' # Generate simulated data for demonstration
#' set.seed(123)
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- rnorm(100)
#'
#' # Fit a Slab Regression model
#' fit <- savvySh(x, y, model_class = "Slab")
#'
#' # Print a detailed summary for all estimators (SR and GSR)
#' summary(fit)
#'
#' # Summarize only a specific estimator
#' summary(fit, estimator = "GSR")
#'
#' @author
#' Ziwei Chen, Vali Asimit, Marina Anca Cidota, Jennifer Asimit\cr
#' Maintainer: Ziwei Chen <ziwei.chen.3@citystgeorges.ac.uk>
#'
#' @importFrom stats pf predict pt qnorm quantile qt
#' @method summary savvySh_model
#' @export
summary.savvySh_model <- function(object, estimator = NULL, ...) {
  cat("Summary of savvySh Model\n")
  cat("===================================================================\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")

  x <- as.matrix(object$model[-1])
  y <- object$model[, 1]

  if (is.null(estimator)) {
    estimator <- names(object$coefficients)
  } else {
    if (!is.character(estimator)) {
      stop("Estimator must be a character vector of estimator names.")
    }
    if (!all(estimator %in% names(object$coefficients))) {
      stop("The specified estimator(s) '", paste(estimator, collapse = ", "),
           "' are not available. Choose from: ", paste(names(object$coefficients), collapse = ", "))
    }
  }

  summary_list <- lapply(estimator, function(est) {
    coefs <- object$coefficients[[est]]
    intercept_present <- (length(coefs) == (ncol(x) + 1))
    stats <- summaryStats_savvySh(x, y, coefs)
    coef_names <- if (intercept_present) {
      c("(Intercept)", colnames(x))
    } else {
      if(is.null(colnames(x))) paste0("V", seq_len(length(coefs))) else colnames(x)
    }
    coef_table <- data.frame(
      Estimate = round(coefs, 4),
      `Std. Error` = round(stats$std_err, 4),
      `t value` = round(stats$t_values, 4),
      `Pr(>|t|)` = format_p_values(stats$p_values),
      `2.5 %` = round(stats$confint_lower, 4),
      `97.5 %` = round(stats$confint_upper, 4),
      `Signif.` = add_significance_codes(stats$p_values),
      check.names = FALSE
    )
    rownames(coef_table) <- coef_names

    list(
      Estimator = est,
      Residuals = stats$residual_quants,
      Coefficients = coef_table,
      ResidualSE = round(stats$residual_se, 4),
      DfResidual = stats$df_residual,
      R_squared = round(stats$r_squared, 4),
      Adj_R_squared = round(stats$adj_r_squared, 4),
      F_statistic = round(stats$f_statistic, 4),
      F_p_value = format_p_values(stats$f_p_value),
      AIC = round(stats$AIC, 4),
      BIC = round(stats$BIC, 4),
      Deviance = round(stats$deviance, 4)
    )
  })

  for (sum_i in summary_list) {
    cat("\nEstimator:", sum_i$Estimator, "\n")
    cat("\nResiduals:\n")
    print(sum_i$Residuals)
    cat("\nCoefficients:\n")
    print(sum_i$Coefficients, row.names = TRUE)
    cat("\nResidual standard error:", sum_i$ResidualSE, "on", sum_i$DfResidual, "degrees of freedom\n")
    cat("Multiple R-squared:", sum_i$R_squared,
        ", Adjusted R-squared:", sum_i$Adj_R_squared, "\n")
    cat("F-statistic:", sum_i$F_statistic,
        "on", (ncol(object$model) - 1), "and", sum_i$DfResidual, "DF, p-value:", sum_i$F_p_value, "\n")
    cat("AIC:", sum_i$AIC, ", BIC:", sum_i$BIC, ", Deviance:", sum_i$Deviance, "\n")
    cat("===================================================================\n")
  }
}
