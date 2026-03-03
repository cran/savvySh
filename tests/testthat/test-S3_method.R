library(testthat)
library(savvySh)

set.seed(123)
x <- matrix(rnorm(100 * 10), ncol = 10)
y <- rnorm(100)
fit <- savvySh(x, y, model_class = "Multiplicative")
fit2 <- savvySh(x, y, model_class = "Linear")

test_that("Print method works correctly for savvySh_model", {
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "Model Class:  Multiplicative")
  expect_output(print(fit), "Estimator:  St")
  expect_output(print(fit), "Estimator:  DSh")
  expect_output(print(fit), "Coefficients")
})

test_that("Summary method works correctly for savvySh_model", {
  expect_output(summary(fit), "Summary of savvySh Model")
  expect_output(summary(fit), "Residuals")
  expect_output(summary(fit), "Coefficients")
  expect_output(summary(fit), "Residual standard error")
  expect_output(summary(fit), "Multiple R-squared")
  expect_output(summary(fit), "F-statistic")
  expect_output(summary(fit), "AIC")
  expect_output(summary(fit), "BIC")
  expect_output(summary(fit), "Deviance")
})

test_that("Coef method works correctly for savvySh_model", {
  coef_St <- coef(fit, estimator = "St")
  coef_DSh <- coef(fit, estimator = "DSh")
  expect_type(coef_St, "double")
  expect_type(coef_DSh, "double")
  expect_warning(coef(fit), "No estimator specified. Returning results for all available estimators: St, DSh")
  expect_equal(length(coef_St), ncol(x) + 1)
  expect_equal(length(coef_DSh), ncol(x) + 1)
})

test_that("Predict method works correctly for savvySh_model", {
  predictions_St <- predict(fit, newx = x[1:5, ], type = "response", estimator = "St")
  predictions_DSh <- predict(fit, newx = x[1:5, ], type = "response", estimator = "DSh")
  predictions_LSh <- predict(fit2, newx = x[1:5, ], type = "response")
  expect_type(predictions_St, "double")
  expect_type(predictions_DSh, "double")
  expect_type(predictions_LSh, "double")
  expect_equal(length(predictions_St), 5)
  expect_equal(length(predictions_DSh), 5)
  expect_equal(length(predictions_LSh), 5)
  coef_St <- predict(fit, type = "coefficients", estimator = "St")
  coef_DSh <- predict(fit, type = "coefficients", estimator = "DSh")
  coef_LSh <- predict(fit2, type = "coefficients")
  expect_type(coef_St, "double")
  expect_type(coef_DSh, "double")
  expect_type(coef_LSh, "double")
  expect_equal(length(coef_St), ncol(x) + 1)
  expect_equal(length(coef_DSh), ncol(x) + 1)
  expect_equal(length(coef_LSh), ncol(x))
  fit_no_intercept <- fit
  fit_no_intercept$coefficients$St <- fit$coefficients$St[-1]
  coef_no_intercept <- predict(fit_no_intercept, type = "coefficients", estimator = "St")
  expect_type(coef_no_intercept, "double")
  expect_equal(length(coef_no_intercept), ncol(x))
  expect_named(coef_no_intercept, paste0("V", seq_len(ncol(x))))
})

test_that("Print method handles invalid and non-standard savvySh_model objects correctly", {
  invalid_fit <- list(
    call = quote(savvySh(x = x, y = y, model_class = "Multiplicative")),
    model_class = "Multiplicative",
    coefficients = NULL
  )
  class(invalid_fit) <- "savvySh_model"
  expect_error(
    capture.output(print(invalid_fit)),
    "The 'savvySh_model' object does not contain valid coefficients."
  )
  invalid_fit$coefficients <- 1
  expect_error(
    capture.output(print(invalid_fit)),
    "The 'savvySh_model' object does not contain valid coefficients."
  )
  expect_error(
    capture.output(print(fit2, estimator = 123)),
    "Estimator must be a character vector"
  )
  expect_error(
    capture.output(print(fit2, estimator = "InvalidEstimator")),
    "The specified estimator\\(s\\) 'InvalidEstimator' are not available. Choose from: LSh"
  )
  no_intercept_fit <- list(
    call = quote(savvySh(x = x, y = y, model_class = "Multiplicative")),
    model_class = "Multiplicative",
    coefficients = list(St = c(0.5, 1.5, 2.5)),
    model = x
  )
  class(no_intercept_fit) <- "savvySh_model"
  output_no_intercept <- capture.output(print(no_intercept_fit))
  expect_true(any(grepl("V1", output_no_intercept)))
  expect_true(any(grepl("V2", output_no_intercept)))
  expect_true(any(grepl("V3", output_no_intercept)))
  expect_false(any(grepl("Intercept", output_no_intercept)))
})

test_that("summary errors when estimator is non-standard", {
  expect_error(
    capture.output(summary(fit2, estimator = 123)),
    "Estimator must be a character vector"
  )
  expect_error(
    capture.output(summary(fit2, estimator = "InvalidEstimator")),
    "The specified estimator\\(s\\) 'InvalidEstimator' are not available. Choose from: LSh"
  )
})

test_that("summary prints coefficients correctly for models without an intercept", {
  out2 <- capture.output(summary(fit2))
  expect_false(any(grepl("\\(Intercept\\)", out2)))
  expect_true(any(grepl("V1", out2)))
})

test_that("Predict method handles missing or improperly structured coefficients", {
  fit_invalid <- fit
  fit_invalid$coefficients <- NULL
  expect_error(
    predict(fit_invalid, newx = x[1:5, ], type = "response", estimator = "St"),
    "Invalid 'savvySh_model' object: coefficients are missing or improperly structured."
  )
  fit_invalid$coefficients <- 1:5
  expect_error(
    predict(fit_invalid, newx = x[1:5, ], type = "response", estimator = "St"),
    "Invalid 'savvySh_model' object: coefficients are missing or improperly structured."
  )
  predictions_default <- suppressWarnings(
    predict(fit, newx = x[1:5, ], type = "response")
  )
  expect_warning(
    predict(fit, newx = x[1:5, ], type = "response"),
    "No estimator specified. Returning results for all available estimators: St, DSh"
  )
  expect_type(predictions_default, "list")
  expect_equal(length(predictions_default$St), 5)
  expect_error(
    predict(fit, newx = x[1:5, ], type = "response", estimator = "InvalidEstimator"),
    "The specified estimator\\(s\\) 'InvalidEstimator' are not available. Choose from: St, DSh"
  )
  fit_invalid$coefficients <- list(St = NULL, DSh = c(0.5, 1.5))
  expect_error(
    predict(fit_invalid, type = "coefficients", estimator = "St"),
    "The coefficients for the estimator St are missing or empty."
  )
  fit_invalid$coefficients <- list(St = numeric(0), DSh = c(0.5, 1.5))
  expect_error(
    predict(fit_invalid, type = "coefficients", estimator = "St"),
    "The coefficients for the estimator St are missing or empty."
  )
})

test_that("Predict method handles errors for 'response' type", {
  expect_error(
    predict(fit, newx = NULL, type = "response", estimator = "St"),
    "You need to supply a matrix of new values for 'newx'."
  )
  newx_invalid <- x[1:5, 1:(ncol(x) - 1)]
  expect_error(
    predict(fit, newx = newx_invalid, type = "response", estimator = "St"),
    "The number of columns in 'newx' must match the number of predictors in the model."
  )
  expect_error(
    predict(fit, newx = x[1:5, ], type = "invalidType", estimator = "St"),
    "Invalid type specified. Use 'response' or 'coefficients'."
  )
})

test_that("predict errors if type is not of length 1 or is not a character vector", {
  expect_error(
    predict(fit, newx = x[1:5, ], type = c("response", "coefficients"), estimator = "St"),
    "Please specify exactly one type: 'response' or 'coefficients'."
  )
  expect_error(
    predict(fit, newx = x[1:5, ], type = "response", estimator = 123),
    "Estimator must be a character vector of estimator names."
  )
})

