library(testthat)
library(savvySh)

set.seed(123)
x <- matrix(rnorm(100 * 10), ncol = 10)
y <- rnorm(100)

test_that("savvySh function handles input validation correctly", {
  expect_error(savvySh(x, y[1:50], model_class = "Multiplicative"), "The number of rows in x must match the length of y.")
  x_with_na <- x
  x_with_na[1, 1] <- NA
  expect_error(savvySh(x_with_na, y, model_class = "Multiplicative"), "x or y contains missing values. Please impute or handle missing data before proceeding.")

  y_with_na <- y
  y_with_na[1] <- NA
  expect_error(savvySh(x, y_with_na, model_class = "Multiplicative"), "x or y contains missing values. Please impute or handle missing data before proceeding.")
  warnings_out <- capture_warnings(
    savvySh(x[1:5, ], y[1:5], model_class = "Multiplicative")
  )

  expect_match(warnings_out[1], "Number of features in x is less than the number of observations.")
  expect_match(warnings_out[2], "Multicollinearity detected: switched to RR instead of unbiased OLS estimation.")
  expect_match(warnings_out[3], "Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold")
})

test_that("savvySh function handles exclusion parameter validation correctly", {
  expect_error(savvySh(x, y, exclude = c(0, 1), model_class = "Multiplicative"), "Exclusion must not contain NA or zero values.")
  expect_error(savvySh(x, y, exclude = NA, model_class = "Multiplicative"), "Exclusion must not contain NA or zero values.")
  expect_error(savvySh(x, y, exclude = c(100), model_class = "Multiplicative"), "Exclusion indices are out of bounds.")
  expect_error(savvySh(x, y, exclude = list(1, 2), model_class = "Multiplicative"), "Exclusion must be a vector.")
  expect_error(savvySh(x, y, exclude = "non_existent_col", model_class = "Multiplicative"), "Exclusion names must match column names.")
  expect_silent(savvySh(x, y, exclude = c(1, 2), model_class = "Multiplicative"))
  expect_silent(savvySh(x, y, exclude = c("V1", "V2"), model_class = "Multiplicative"))
  expect_error(savvySh(x, y, exclude = c(1, "V2"), model_class = "Multiplicative"), "Exclusion names must match column names.")
})

test_that("savvySh function validates model_class argument correctly", {
  valid_models <- c("Multiplicative", "Slab", "Linear", "ShrinkageRR")
  for (model in valid_models) {
    expect_silent(savvySh(x, y, model_class = model))
  }
  expect_error(
    savvySh(x, y, model_class = "InvalidModel"),
    "Invalid model_class specified. Choose from: Multiplicative, Slab, Linear, ShrinkageRR"
  )
  expect_warning(savvySh(x, y, model_class = "Slab", include_Sh = TRUE),
                 "include_Sh is only applicable when model_class is 'Multiplicative'. Setting include_Sh to FALSE.")
})

test_that("savvySh warns and uses the first model when multiple model_class values are provided", {
  warn_msg <- capture_warnings(savvySh(x, y, model_class = c("Multiplicative", "Slab")))
  expect_match(warn_msg, "model_class should be a single option. Using the first element:")
  result <- suppressWarnings(savvySh(x, y, model_class = c("Multiplicative", "Slab")))
  expect_equal(result$model_class, "Multiplicative")
})

test_that("savvySh function validates lambda_vals and folds arguments correctly", {
  x_multicollinear <- cbind(x, x[, 1])
  lambda_vals <- 10^seq(-3, 1, length.out = 5)
  expect_error(
    savvySh(x_multicollinear, y, model_class = "Multiplicative", lambda_vals = c(0.1), folds = 5),
    "Need more than one value of lambda_vals for cross-validation."
  )
  expect_error(
    savvySh(x_multicollinear, y, model_class = "Multiplicative", lambda_vals = lambda_vals, folds = 2.5),
    "Number of folds must be an integer greater than or equal to 3."
  )
  expect_error(
    savvySh(x_multicollinear, y, model_class = "Multiplicative", lambda_vals = lambda_vals, folds = 2),
    "Number of folds must be an integer greater than or equal to 3."
  )
})

test_that("savvySh function works correctly for Multiplicative model with and without Sh estimator", {
  result_mult_no_sh <- savvySh(x, y, model_class = "Multiplicative", include_Sh = FALSE)
  expect_s3_class(result_mult_no_sh, "savvySh_model")
  expect_named(result_mult_no_sh, c("call", "model", "model_class", "optimal_lambda", "coefficients", "fitted_values", "pred_MSE"))
  expect_named(result_mult_no_sh$coefficients, c("St", "DSh"))
  expect_named(result_mult_no_sh$fitted_values, c("St", "DSh"))
  expect_named(result_mult_no_sh$pred_MSE, c("St", "DSh"))

  result_mult_with_sh <- savvySh(x, y, model_class = "Multiplicative", include_Sh = TRUE)
  expect_named(result_mult_with_sh$coefficients, c("St", "DSh", "Sh"))
  expect_named(result_mult_with_sh$fitted_values, c("St", "DSh", "Sh"))
  expect_named(result_mult_with_sh$pred_MSE, c("St", "DSh", "Sh"))
})

test_that("savvySh function works correctly for Slab model with valid v parameter", {
  result_slab <- savvySh(x, y, model_class = "Slab", v = 2)
  expect_s3_class(result_slab, "savvySh_model")
  expect_named(result_slab, c("call", "model", "model_class", "optimal_lambda", "coefficients", "fitted_values", "pred_MSE"))
  expect_named(result_slab$coefficients, c("SR", "GSR"))
  expect_named(result_slab$fitted_values, c("SR", "GSR"))
  expect_named(result_slab$pred_MSE, c("SR", "GSR"))
})

test_that("savvySh function throws error for invalid v parameter in Slab model", {
  expect_error(savvySh(x, y, model_class = "Slab", v = -1), "Error: 'v' must be a positive number.")
  expect_error(savvySh(x, y, model_class = "Slab", v = "non-numeric"), "Error: 'v' must be a positive number.")
})

test_that("savvySh function works correctly for Linear model with automatic centering", {
  result_linear <- savvySh(x, y, model_class = "Linear")
  expect_s3_class(result_linear, "savvySh_model")
  expect_named(result_linear, c("call", "model", "model_class", "optimal_lambda", "coefficients", "fitted_values", "pred_MSE"))
  expect_named(result_linear$coefficients, "LSh")
  expect_named(result_linear$fitted_values, "LSh")
  expect_named(result_linear$pred_MSE, "LSh")

  centered_x <- scale(x, center = TRUE, scale = FALSE)
  centered_y <- scale(y, center = TRUE, scale = FALSE)
  result_centered <- savvySh(centered_x, centered_y, model_class = "Linear")
  expect_equal(result_linear$coefficients$LSh, result_centered$coefficients$LSh, tolerance = 1e-6)
})

test_that("savvySh function works correctly for SRR model", {
  result_slab <- savvySh(x, y, model_class = "ShrinkageRR")
  expect_s3_class(result_slab, "savvySh_model")
  expect_named(result_slab, c("call", "model", "model_class", "coefficients", "fitted_values", "pred_MSE"))
  expect_named(result_slab$coefficients, "SRR")
  expect_named(result_slab$fitted_values, "SRR")
  expect_named(result_slab$pred_MSE, "SRR")
})

test_that("savvySh function requires at least two covariates for Linear model", {
  x_single_covariate <- matrix(rnorm(100), ncol = 1)
  y <- rnorm(100)
  expect_error(
    savvySh(x_single_covariate, y, model_class = "Linear"),
    "The number of covariates must be at least two for this model."
  )
})

test_that("savvySh function handles multicollinearity correctly with RR in Multiplicative model", {
  x_multicollinear <- cbind(x, x[, 1])
  expect_warning(
    result_multicollinear_mult <- savvySh(x_multicollinear, y, model_class = "Multiplicative"),
    "Multicollinearity detected: switched to RR instead of unbiased OLS estimation."
  )
  expect_true("ridge_results" %in% names(result_multicollinear_mult))
  expect_named(result_multicollinear_mult$ridge_results, c("lambda_range", "cvm", "cvsd", "ridge_coefficients"))
  expect_equal(length(result_multicollinear_mult$ridge_results$lambda_range), 100)
  expect_true(is.numeric(result_multicollinear_mult$ridge_results$lambda_range))
  expect_true(is.numeric(result_multicollinear_mult$ridge_results$cvm))
  expect_true(is.numeric(result_multicollinear_mult$ridge_results$cvsd))
  expect_true(is.matrix(result_multicollinear_mult$ridge_results$ridge_coefficients))
})

test_that("savvySh function handles multicollinearity correctly with RR in Linear model", {
  x_multicollinear <- cbind(x, x[, 1])
  lambda_vals <- 10^seq(-3, 2, length.out = 5)
  expect_warning(
    result_multicollinear_linear <- savvySh(x_multicollinear, y, model_class = "Linear", lambda_vals = lambda_vals, folds = 5, foldid = TRUE),
    "Multicollinearity detected: switched to RR on centered data instead of unbiased OLS estimation."
  )
  expect_true("ridge_results" %in% names(result_multicollinear_linear))
  expect_named(result_multicollinear_linear$ridge_results, c("lambda_range", "cvm", "cvsd", "ridge_coefficients"))
  expect_true("fold_assignments" %in% names(result_multicollinear_linear))
  expect_length(result_multicollinear_linear$fold_assignments, nrow(x))
  expect_equal(sort(result_multicollinear_linear$ridge_results$lambda_range), sort(lambda_vals))
  expect_true(is.numeric(result_multicollinear_linear$ridge_results$lambda_range))
  expect_true(is.numeric(result_multicollinear_linear$ridge_results$cvm))
  expect_true(is.numeric(result_multicollinear_linear$ridge_results$cvsd))
  expect_true(is.matrix(result_multicollinear_linear$ridge_results$ridge_coefficients))
})


