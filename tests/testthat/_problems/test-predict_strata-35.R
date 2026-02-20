# Extracted from test-predict_strata.R:35

# test -------------------------------------------------------------------------
set.seed(42)
x <- rlnorm(200, meanlog = 6, sdlog = 1.5)
sb <- strata_bound(x, n_strata = 3, n = 60)
newx <- c(-1, x, max(x) + 1e6)
f <- predict(sb, newx)
expect_false(anyNA(f))
