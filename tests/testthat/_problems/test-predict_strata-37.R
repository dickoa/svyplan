# Extracted from test-predict_strata.R:37

# test -------------------------------------------------------------------------
set.seed(42)
x <- rlnorm(200, meanlog = 6, sdlog = 1.5)
sb <- strata_bound(x, n_strata = 3, n = 60)
newx <- c(-1, x, max(x) + 1e6)
f <- predict(sb, newx)
expect_false(anyNA(f))
expect_equal(as.integer(f[1L]), 1L)
expect_equal(as.integer(f[length(f)]), 3L)
