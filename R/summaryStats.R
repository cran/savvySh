#' @keywords internal
format_p_values <- function(p_values, digits = 4) {
  formatted_p_values <- ifelse(
    p_values < 10^(-digits),
    sprintf(paste0("%.", digits, "e"), p_values),
    sprintf(paste0("%.", digits, "f"), p_values)
  )
  return(formatted_p_values)
}

#' @keywords internal
add_significance_codes <- function(p_values) {
  signif_codes <- rep("", length(p_values))
  signif_codes[p_values <= 0.001] <- "***"
  signif_codes[p_values > 0.001 & p_values <= 0.01] <- "**"
  signif_codes[p_values > 0.01 & p_values <= 0.05] <- "*"
  signif_codes[p_values > 0.05 & p_values <= 0.1] <- "."
  signif_codes[p_values > 0.1] <- " "
  return(signif_codes)
}

#' @keywords internal
#' Compute summary statistics for a given set of coefficients in a linear model.
#' This is used by \code{summary.savvySh_model} to provide OLS-like diagnostics,
#' including residuals, R-squared, F-stats, and confidence intervals.
#'
summaryStats_savvySh <- function(x, y, coefficients, conf_level = 0.95) {
  nobs <- nrow(x)
  nvars <- ncol(x)

  intercept_present <- length(coefficients) == (nvars + 1)
  if (intercept_present) {
    x <- cbind(1, x)
  }

  fitted_values <- x %*% coefficients
  residuals <- y - fitted_values
  residual_quants <- quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1))

  df_residual <- nobs - (nvars + as.integer(intercept_present))
  rss <- sum(residuals^2)
  residual_se <- sqrt(rss / df_residual)

  rss <- sum(residuals^2)
  df_residual <- nobs - nvars - as.integer(intercept_present)
  residual_se <- sqrt(rss / df_residual)

  tss <- sum((y - mean(y))^2)
  r_squared <- 1 - (rss / tss)
  adj_r_squared <- 1 - (1 - r_squared) * ((nobs - 1) / df_residual)

  xtx_inv <- solve(crossprod(x))
  std_err <- residual_se * sqrt(diag(xtx_inv))

  t_values <- coefficients / std_err
  p_values <- 2 * pt(-abs(t_values), df = df_residual)

  msr <- (tss - rss) / nvars
  mse <- rss / df_residual
  f_statistic <- msr / mse
  f_p_value <- pf(f_statistic, nvars, df_residual, lower.tail = FALSE)

  log_likelihood <- -0.5 * nobs * log(rss / nobs)
  k <- nvars + as.integer(intercept_present)
  AIC <- -2 * log_likelihood + 2 * k
  BIC <- -2 * log_likelihood + log(nobs) * k
  deviance <- rss

  alpha <- 1 - conf_level
  t_val <- qt(1 - alpha / 2, df = df_residual)
  confint_lower <- coefficients - t_val * std_err
  confint_upper <- coefficients + t_val * std_err

  list(
    residuals = residuals,
    residual_quants = residual_quants,
    residual_se = residual_se,
    std_err = std_err,
    t_values = t_values,
    p_values = p_values,
    df_residual = df_residual,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    f_statistic = f_statistic,
    f_p_value = f_p_value,
    AIC = AIC,
    BIC = BIC,
    deviance = deviance,
    confint_lower = confint_lower,
    confint_upper = confint_upper
  )
}

