# savvySh: Shrinkage Methods for Linear Regression Estimation

The `savvySh` package provides a unified interface for fitting shrinkage estimators in linear regression, which is particularly useful in the presence of multicollinearity or high-dimensional covariates. It supports four shrinkage classes: Multiplicative Shrinkage, Slab Regression, Linear Shrinkage, and Shrinkage Ridge Regression. These methods improve on the classical Ordinary Least Squares (OLS) estimator by trading a small amount of bias for a significant reduction in variance.

This package implements the theoretical framework discussed in:

Asimit, V., Cidota, M. A., Chen, Z., & Asimit, J. (2025). [Slab and Shrinkage Linear Regression Estimation](https://openaccess.city.ac.uk/id/eprint/35005/).

## Related Projects

-   **savvyGLM**: For applying these shrinkage methods within Generalized Linear Models (GLMs), please refer to the companion package [savvyGLM](https://github.com/Ziwei-ChenChen/savvyGLM).
-   **flashfm-savvySh**: For applications in genetic fine-mapping, see the [flashfm-savvySh](https://github.com/jennasimit/flashfm-savvySh) repository.

## Installation Guide

You can install the released version of `savvySh` from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("savvySh")
```

Alternatively, you can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("Ziwei-ChenChen/savvySh")
```

Once installed, load the package:

``` r
library(savvySh)
```

## Features

`savvySh` provides several shrinkage estimators designed to improve regression accuracy by reducing Mean Squared Error (MSE):

-   **Multiplicative Shrinkage:** Applies shrinkage by multiplying the OLS estimates with data-driven factors.

    -   **Stein (St):** Applies a single global shrinkage factor to all coefficients.

    -   **Diagonal Shrinkage (DSh):** Applies a separate factor to each coefficient.

    -   **Shrinkage (Sh):** Uses a full matrix shrinkage operator estimated by solving a *Sylvester equation*.

-   **Slab Regression:** Adds structured shrinkage based on penalty terms.

    -   **Slab Regression (SR):** Shrinks toward a fixed target direction (e.g., a vector of ones).

    -   **Generalized Slab Regression (GSR):** Shrinks toward multiple directions (e.g., eigenvectors).

-   **Linear Shrinkage (LSh):** Takes a weighted average of the OLS estimator and a target estimator and is useful for standardized data.

-   **Shrinkage Ridge Regression (SRR):** Extends *Ridge Regression (RR)* by shrinking toward a diagonal matrix with equal entries.

All shrinkage factors are computed in closed form (except SRR, which optimizes shrinkage intensity numerically).

## Usage

This is a basic example that shows you how to solve a common problem:

``` r
# Simulated example
set.seed(123)
x <- matrix(rnorm(100 * 10), 100, 10)
y <- rnorm(100)

# Fit shrinkage estimators
fit <- savvySh(x, y, model_class = "Multiplicative", include_Sh = TRUE)

# Extract coefficients
coef(fit, estimator = "St")
coef(fit, estimator = "DSh")
coef(fit, estimator = "Sh")
```

## Authors

-   Ziwei Chen – [ziwei.chen.3\@citystgeorges.ac.uk](mailto:ziwei.chen.3@citystgeorges.ac.uk)
-   Vali Asimit – [asimit\@citystgeorges.ac.uk](mailto:asimit@citystgeorges.ac.uk)
-   Marina Anca Cidota – [cidota\@fmi.unibuc.ro](mailto:cidota@fmi.unibuc.ro)
-   Jennifer Asimit – [jennifer.asimit\@mrc-bsu.cam.ac.uk](mailto:jennifer.asimit@mrc-bsu.cam.ac.uk)

## License

This package is licensed under the GPL (\>= 3) License.
