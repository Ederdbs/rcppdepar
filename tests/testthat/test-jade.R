## JADE (control$c > 0) adaptation of CR / F.
##
## devol() used to declare goodCR / goodF / goodF2 / i_goodNP outside the
## generation loop and never reset them, so from the second generation on the
## adaptation was driven by all-time accumulations. Worse, it computed
## goodF2/goodF unguarded: a generation with no successful trial made that
## 0.0/0.0 = NaN, meanF went NaN, R::rcauchy(NaN, 0.1) returned NaN, the
## do/while accepted it (NaN comparisons are false), and every trial vector
## after that was NaN -- which the bounce-back bounds check does not catch
## either, since NaN < lower and NaN > upper are both false.

test_that("a generation with zero successful trials does not poison meanF", {
  ## This objective makes acceptance impossible after initialization: the
  ## first NP calls (the initial population) score 0, every mutant scores
  ## sum(x^2) >= 0, and acceptance needs strictly-better-or-equal. So nGood
  ## is 0 in *every* generation -- the exact 0/0 case.
  ##
  ## Under the old code the NaN weight produced NaN trial vectors, and
  ## EvalStandard raised "NaN value of objective function!".
  NP <- 10
  n <- 0
  fn <- function(x) {
    n <<- n + 1
    if (n <= NP) 0 else sum(x^2)
  }

  set.seed(11)
  res <- suppressWarnings(
    DEoptim(fn, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = NP, itermax = 25, trace = FALSE, c = 0.5)))

  expect_true(is.finite(res$optim$bestval))
  expect_equal(res$optim$bestval, 0)
  expect_true(all(is.finite(res$member$pop)))
  expect_true(all(is.finite(res$member$bestvalit)))
})

test_that("JADE runs stay finite across all strategies", {
  for (s in 1:6) {
    set.seed(100 + s)
    res <- suppressWarnings(
      DEoptim(function(x) sum(x^2), rep(-5, 3), rep(5, 3),
              DEoptim.control(NP = 20, itermax = 40, trace = FALSE,
                              strategy = s, c = 0.5)))
    expect_true(is.finite(res$optim$bestval), info = paste("strategy", s))
    expect_true(all(is.finite(res$member$bestvalit)), info = paste("strategy", s))
    expect_true(all(is.finite(res$member$pop)), info = paste("strategy", s))
  }
})

test_that("JADE still optimizes (adaptation is not frozen or NaN)", {
  set.seed(5)
  res <- suppressWarnings(
    DEoptim(function(x) sum(x^2), rep(-10, 5), rep(10, 5),
            DEoptim.control(NP = 50, itermax = 200, trace = FALSE, c = 0.5)))
  expect_lt(res$optim$bestval, 1e-6)
  ## bestvalit is monotone non-increasing by construction; a NaN anywhere in
  ## the adaptation would break that.
  expect_true(all(diff(res$member$bestvalit) <= 0))
})

test_that("c == 0 disables JADE and still converges", {
  set.seed(5)
  res <- suppressWarnings(
    DEoptim(function(x) sum(x^2), rep(-10, 5), rep(10, 5),
            DEoptim.control(NP = 50, itermax = 200, trace = FALSE, c = 0)))
  expect_lt(res$optim$bestval, 1e-6)
})
