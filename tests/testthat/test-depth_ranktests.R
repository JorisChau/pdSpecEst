context("Geometric functions")

test_that("Correctly working data depth and rank-based tests", {

  # Data depth
  d <- 2
  mu <- diag(d)
  imu <- replicate(3, mu)
  E <- pdSpecEst:::E_basis(d)

  ## Pointwise depth
  X1 <- replicate(10, Expm(mu, pdSpecEst:::E_coeff_inv(rnorm(d^2), E)))
  dd1 <- pdDepth(mu, X1, "zonoid")
  dd2 <- pdDepth(mu, X1, "gdd")
  dd3 <- pdDepth(mu, X1, "spatial")
  dd4 <- pdDepth(X = X1, method = "zonoid")
  dd5 <- pdDepth(X = X1, method = "gdd")
  dd6 <- pdDepth(X = X1, method = "spatial")

  ## Integrated depth
  X2 <- replicate(10, replicate(3, Expm(mu, pdSpecEst:::E_coeff_inv(rnorm(d^2), E))))
  idd1 <- pdDepth(imu, X2, "zonoid")
  idd2 <- pdDepth(imu, X2, "gdd")
  idd3 <- pdDepth(imu, X2, "spatial")
  idd4 <- pdDepth(X = X2, method = "zonoid")
  idd5 <- pdDepth(X = X2, method = "gdd")
  idd6 <- pdDepth(X = X2, method = "spatial")

  expect_true(is.numeric(dd1) & (length(dd1) == 1))
  expect_true(is.numeric(dd2) & (length(dd2) == 1))
  expect_true(is.numeric(dd3) & (length(dd3) == 1))
  expect_true(is.numeric(idd1) & (length(idd1) == 1))
  expect_true(is.numeric(idd2) & (length(idd2) == 1))
  expect_true(is.numeric(idd3) & (length(idd3) == 1))

  expect_true(is.numeric(dd4) & (length(dd4) == 10))
  expect_true(is.numeric(dd5) & (length(dd5) == 10))
  expect_true(is.numeric(dd6) & (length(dd6) == 10))
  expect_true(is.numeric(idd4) & (length(idd4) == 10))
  expect_true(is.numeric(idd5) & (length(idd5) == 10))
  expect_true(is.numeric(idd6) & (length(idd6) == 10))

  # Rank-based tests
  rs <- pdRankTests(X1, c(5, 5), "rank.sum")
  irs <- pdRankTests(X2, c(5, 5), "rank.sum")
  kw <- pdRankTests(X1, c(5, 3, 2), "krusk.wall")
  ikw <- pdRankTests(X2, c(5, 3, 2), "krusk.wall")
  sr <- pdRankTests(X1, test = "signed.rank")
  bvn <- pdRankTests(X1, test = "bartels")
  ibvn <- pdRankTests(X2, test = "bartels")

  expect_true(is.list(rs) & (length(rs) == 3))
  expect_match(rs[[3]], "Standard normal distribution")
  expect_true(is.list(irs) & (length(irs) == 3))
  expect_match(irs[[3]], "Standard normal distribution")
  expect_true(is.list(kw) & (length(kw) == 3))
  expect_true(is.list(ikw) & (length(ikw) == 3))
  expect_true(is.list(sr) & (length(sr) == 3))
  expect_match(sr[[3]], "Wilcoxon signed rank test")
  expect_true(is.list(bvn) & (length(bvn) == 3))
  expect_match(bvn[[3]], "Standard normal distribution")
  expect_true(is.list(ibvn) & (length(ibvn) == 3))
  expect_match(ibvn[[3]], "Standard normal distribution")

})
