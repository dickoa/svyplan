test_that("range moments agree with direct statistics across scales", {
  populations <- list(
    1:100, 1e9 + 1:100, 1e-12*(1:100),
    c(1:50, 1e12 + 1:50),
    c(1:50, rep(1e12,50)),
    c(rep(1,50),1e12 + (1:50)/100),
    rep(1e9,100)
  )
  for (x in populations) {
    idx <- c(0,1,20,50,80,99,100)
    st <- .strata_stats_from_prefix(.strata_precompute(x),idx)
    for (j in seq_len(length(idx)-1L)) {
      values <- x[seq.int(idx[j]+1,idx[j+1])]
      expect_equal(st$mean_h[j],mean(values),tolerance=1e-12)
      want <- if(length(values)==1L) 0 else sd(values)
      expect_equal(st$S_h[j],want,tolerance=1e-8)
    }
  }
  # A singleton has no variance to trigger the cancellation fallback, and
  # reconstructing its tiny mean from the global center can lose it entirely.
  x <- c(1e-15, 1e12 + 1:99)
  st <- .strata_stats_from_prefix(.strata_precompute(x), c(0,1,100))
  expect_identical(st$mean_h[1L],x[1L])
})

test_that("a large location does not change the fixed-n LH design", {
  a <- strata_bound(1:100,n_strata=2,n=20,method="lh")
  b <- strata_bound(1e9+1:100,n_strata=2,n=20,method="lh")
  expect_equal(a$strata$N,c(50,50))
  expect_equal(b$strata$N,a$strata$N)
  expect_equal(b$strata$n,a$strata$n)
  expect_equal((b$cv*mean(1e9+1:100))^2,8.5,tolerance=1e-10)
})

test_that("both searches preserve designs under translation", {
  x <- (1:120)^2
  for (method in c("lh","kozak")) {
    set.seed(907)
    a <- strata_bound(x,n_strata=3,n=30,method=method)
    set.seed(907)
    b <- strata_bound(1e9+x,n_strata=3,n=30,method=method)
    expect_equal(b$strata$N,a$strata$N,info=method)
    expect_equal(b$strata$n,a$strata$n,info=method)
    expect_equal(b$cv*mean(1e9+x),a$cv*mean(x),tolerance=1e-8)
    set.seed(908)
    c <- strata_bound(x,n_strata=3,cv=.03,method=method)
    set.seed(908)
    d <- strata_bound(1e9+x,n_strata=3,
                      cv=.03*mean(x)/mean(1e9+x),method=method)
    expect_equal(d$strata$N,c$strata$N,info=method)
    expect_equal(d$strata$n,c$strata$n,info=method)
  }
})
