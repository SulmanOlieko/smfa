library(smfa)

test_that("smfa works with sfacross and different metaMethods", {
  # Generate small toy data
  set.seed(123)
  N <- 80
  z1 <- rnorm(N)
  v1 <- rnorm(N)
  g <- rnorm(N)
  ds <- z1 + v1
  d <- ifelse(ds > 0, 1, 0)
  group <- ifelse(g > 0, 1, 0)
  x1 <- rnorm(N)
  y <- x1 + rnorm(N) - abs(rnorm(N))
  dat <- data.frame(y = y, x1 = x1, z1 = z1, d = d, group = group)
  
  # 1. LP Metafrontier
  meta_lp <- smfa(
    formula    = y ~ x1,
    data       = dat,
    group      = "group",
    groupType  = "sfacross",
    metaMethod = "lp"
  )
  expect_s3_class(meta_lp, "smfa")
  expect_equal(meta_lp$metaMethod, "lp")
  
  # 2. QP Metafrontier
  suppressWarnings({
    meta_qp <- smfa(
      formula    = y ~ x1,
      data       = dat,
      group      = "group",
      groupType  = "sfacross",
      metaMethod = "qp"
    )
  })
  expect_s3_class(meta_qp, "smfa")
  expect_equal(meta_qp$metaMethod, "qp")
  
  # 3. SFA Metafrontier - Huang
  suppressWarnings({
    meta_huang <- smfa(
      formula     = y ~ x1,
      data        = dat,
      group       = "group",
      groupType   = "sfacross",
      metaMethod  = "sfa",
      sfaApproach = "huang"
    )
  })
  expect_s3_class(meta_huang, "smfa")
  expect_equal(meta_huang$sfaApproach, "huang")
  
  # Test S3 methods on meta_lp and meta_huang
  expect_output(print(meta_lp))
  expect_output(print(summary(meta_lp)))
  expect_type(coef(meta_lp), "list") # coef returns a list for LP
  expect_s3_class(efficiencies(meta_lp), "data.frame")
  expect_true(is.numeric(fitted(meta_lp)))
  expect_true(is.numeric(residuals(meta_lp)))
  
  ll <- logLik(meta_lp)
  expect_true(is.numeric(as.numeric(ll)))
  expect_true(is.numeric(nobs(meta_lp)))
  expect_s3_class(ic(meta_lp), "data.frame")
  
  # Test coef and ic on a method that returns them (sfa)
  expect_true(is.numeric(coef(meta_huang)))
  expect_s3_class(ic(meta_huang), "data.frame")
})

test_that("smfa works with sfaselectioncross", {
  set.seed(12345)
  N <- 300 # increased N to ensure group-specific sample sizes > 20
  z1 <- rnorm(N)
  v1 <- rnorm(N)
  g <- rnorm(N)
  ds <- z1 + v1
  d <- ifelse(ds > 0, 1, 0)
  group <- ifelse(g > 0, 1, 0)
  x1 <- rnorm(N)
  y <- x1 + rnorm(N) - abs(rnorm(N))
  dat <- data.frame(y = y, x1 = x1, z1 = z1, d = d, group = group)
  
  meta_sel <- smfa(
    formula    = y ~ x1,
    selectionF = d ~ z1,
    data       = dat,
    group      = "group",
    groupType  = "sfaselectioncross",
    lType      = "ghermite",
    Nsub       = 5,
    itermax    = 20,
    metaMethod = "lp"
  )
  expect_s3_class(meta_sel, "smfa")
})
