test_that("Test linear regression.", {
  
  # Validate linear regression.
  withr::local_seed(1010)
  n <- 1000
  x <- rnorm(n = n)
  y <- rnorm(n = n)
  
  fit_lm <- lm(y ~ x)
  fit_pkg <- FitOLS(y, cbind(1, x))
  
  # Compare with built-in functions.
  beta_lm <- coef(fit_lm)
  beta_pkg <- as.numeric(fit_pkg$beta)
  expect_equal(beta_pkg, beta_lm, ignore_attr = TRUE)
  
  summary_lm <- summary(fit_lm)
  se_lm <- summary_lm$coefficients[, "Std. Error"]
  se_pkg <- as.numeric(fit_pkg$se)
  expect_equal(se_pkg, se_lm, ignore_attr = TRUE)
  
  sigma2_lm <- sigma(fit_lm)^2
  sigma2_pkg <- fit_pkg$sigma
  expect_equal(sigma2_pkg, sigma2_lm, ignore_attr = TRUE)
  
})


test_that("Test weighted linear regression.", {
  
  # Validate linear regression.
  withr::local_seed(1010)
  n <- 1000
  x <- rnorm(n = n)
  y <- rnorm(n = n)
  w <- runif(n)
  
  fit_lm <- lm(y ~ x, weights = w)
  fit_pkg <- FitWLS(y, cbind(1, x), w = w)
  
  # Compare with built-in functions.
  beta_lm <- coef(fit_lm)
  beta_pkg <- as.numeric(fit_pkg$beta)
  expect_equal(beta_pkg, beta_lm, ignore_attr = TRUE)
  
  summary_lm <- summary(fit_lm)
  se_lm <- summary_lm$coefficients[, "Std. Error"]
  se_pkg <- as.numeric(fit_pkg$se)
  expect_equal(se_pkg, se_lm, ignore_attr = TRUE)
  
  sigma2_lm <- sigma(fit_lm)^2
  sigma2_pkg <- fit_pkg$sigma
  expect_equal(sigma2_pkg, sigma2_lm, ignore_attr = TRUE)
  
})
