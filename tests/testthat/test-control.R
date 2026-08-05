## Input validation in DEoptim.control() / DEoptim().
## Each check here corresponds to a value that used to reach C++ and either
## overrun an array or divide by zero.

test_that("NP below 6 is rejected (permute() draws 5 distinct indices)", {
  for (np in c(1, 3, 4, 5)) {
    expect_warning(ctrl <- DEoptim.control(NP = np, trace = FALSE), "'NP' < 6")
    expect_equal(ctrl$NP, 50)
  }
  expect_silent(ctrl <- DEoptim.control(NP = 6, trace = FALSE))
  expect_equal(ctrl$NP, 6)
})

test_that("NP == 6 runs without overrunning the urn arrays", {
  ## NP == 4 used to write ia_urn2[5] (array of 5) and ia_urn1[-1] here.
  ## 6 is the new minimum, i.e. the tightest case permute() has to handle.
  set.seed(42)
  res <- suppressWarnings(
    DEoptim(function(x) sum(x^2), rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 6, itermax = 30, trace = FALSE)))
  expect_true(is.finite(res$optim$bestval))
  expect_equal(dim(res$member$pop), c(6L, 2L))
})

test_that("C++ entry point refuses NP < 6 on its own", {
  ## DEoptim.control() clamps, so reach DEoptim_impl() directly with a
  ## hand-built control list.
  ctrl <- suppressWarnings(DEoptim.control(NP = 6, itermax = 5, trace = FALSE))
  ctrl$NP <- 4
  ctrl$npar <- 2
  ctrl$specinitialpop <- 0
  ctrl$initialpop <- matrix(0, 1, 1)
  ctrl$trace <- 0
  expect_error(
    RcppDEpar:::DEoptim_impl(rep(-5, 2), rep(5, 2), function(x) sum(x^2),
                             ctrl, new.env()),
    "'NP' must be at least 6")
})

test_that("storepopfreq below 1 is rejected (used to divide by zero in C++)", {
  expect_warning(ctrl <- DEoptim.control(storepopfreq = 0, trace = FALSE),
                 "'storepopfreq' < 1")
  expect_equal(ctrl$storepopfreq, 1)

  expect_warning(ctrl <- DEoptim.control(storepopfreq = -3, trace = FALSE),
                 "'storepopfreq' < 1")
  expect_equal(ctrl$storepopfreq, 1)
})

test_that("nthreads > 1 with a non-thread-safe objective warns", {
  ## The parallel path is only taken for a raw-tagged compiled objective;
  ## an R function silently ran serially before.
  set.seed(1)
  expect_warning(
    DEoptim(function(x) sum(x^2), rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 20, itermax = 5, trace = FALSE, nthreads = 4)),
    "requires a thread-safe compiled objective")
})
