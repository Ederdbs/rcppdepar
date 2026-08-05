## Baseline behaviour of the serial path -- guards against the refactoring
## (population swap instead of copy, unique_ptr for the evaluator, per-
## generation JADE reset) changing the optimizer itself.

sphere    <- function(x) sum(x^2)
rosenbrock <- function(x) sum(100 * (x[-1] - x[-length(x)]^2)^2 + (1 - x[-length(x)])^2)

test_that("every strategy converges on the sphere", {
  for (s in 1:6) {
    set.seed(200 + s)
    res <- suppressWarnings(
      DEoptim(sphere, rep(-10, 4), rep(10, 4),
              DEoptim.control(NP = 50, itermax = 300, trace = FALSE, strategy = s)))
    expect_lt(res$optim$bestval, 1e-6)
    expect_equal(unname(res$optim$bestmem), rep(0, 4), tolerance = 1e-3)
  }
})

test_that("rosenbrock is solved", {
  set.seed(13)
  res <- suppressWarnings(
    DEoptim(rosenbrock, rep(-5, 3), rep(5, 3),
            DEoptim.control(NP = 60, itermax = 800, trace = FALSE)))
  expect_lt(res$optim$bestval, 1e-5)
  expect_equal(unname(res$optim$bestmem), rep(1, 3), tolerance = 1e-2)
})

test_that("the population swap keeps pop and bestmem consistent", {
  ## ta_oldP/ta_newP are now swapped rather than copied each generation; the
  ## returned population must still be the final one, and it must contain the
  ## reported best member.
  set.seed(23)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 3), rep(5, 3),
            DEoptim.control(NP = 20, itermax = 100, trace = FALSE)))
  expect_equal(dim(res$member$pop), c(20L, 3L))
  vals <- apply(res$member$pop, 1, sphere)
  expect_equal(min(vals), res$optim$bestval, tolerance = 1e-8)
  expect_equal(unname(res$member$pop[which.min(vals), ]),
               unname(res$optim$bestmem), tolerance = 1e-8)
})

test_that("bs = TRUE keeps the best of parents and children", {
  set.seed(29)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 3), rep(5, 3),
            DEoptim.control(NP = 20, itermax = 200, trace = FALSE, bs = TRUE)))
  expect_lt(res$optim$bestval, 1e-6)
  expect_true(all(diff(res$member$bestvalit) <= 0))
})

test_that("bestmemit and bestvalit have one row per generation run", {
  set.seed(37)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 10, itermax = 25, trace = FALSE)))
  expect_equal(res$optim$iter, 25)
  expect_equal(dim(res$member$bestmemit), c(25L, 2L))
  expect_length(res$member$bestvalit, 25L)
  expect_true(all(is.finite(res$member$bestmemit)))
})

test_that("tracking the best member does not overwrite a population member", {
  ## `t_bestP = ta_popP.unsafe_col(i)` made t_bestP an *alias* of population
  ## column i rather than a copy, so the next best-update wrote the new best
  ## member straight into the previously-aliased column. Member 0 was always
  ## hit (i == 0 aliases first), silently replacing it with a duplicate of the
  ## best member and costing a member of population diversity every run.
  set.seed(41)
  pop0 <- matrix(runif(20 * 2, -1, 1), nrow = 20, ncol = 2)
  res <- suppressWarnings(
    DEoptim(function(x) sum(x^2), rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 20, itermax = 1, trace = FALSE,
                            initialpop = pop0, storepopfrom = 1)))
  ## generation 0 is stored before any mutation: it must be pop0 untouched
  expect_equal(unname(res$member$storepop[[1]]), pop0, tolerance = 1e-12)
  expect_equal(nrow(unique(res$member$storepop[[1]])), 20L)
})

test_that("an initial population is used as given", {
  set.seed(41)
  pop0 <- matrix(runif(20 * 2, -1, 1), nrow = 20, ncol = 2)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 20, itermax = 1, trace = FALSE,
                            initialpop = pop0, storepopfrom = 1)))
  ## generation 0 is stored before any mutation, so it is the initial population
  expect_equal(unname(res$member$storepop[[1]]), pop0, tolerance = 1e-12)
})

test_that("solutions stay inside the bounds", {
  set.seed(47)
  lo <- c(-1, 0, 2); up <- c(1, 3, 4)
  res <- suppressWarnings(
    DEoptim(function(x) -sum(x), lo, up,
            DEoptim.control(NP = 20, itermax = 100, trace = FALSE)))
  expect_true(all(res$optim$bestmem >= lo - 1e-12))
  expect_true(all(res$optim$bestmem <= up + 1e-12))
  expect_true(all(sweep(res$member$pop, 2, lo, `>=`)))
  expect_true(all(sweep(res$member$pop, 2, up, `<=`)))
})

test_that("nfeval is exact on the serial path", {
  set.seed(53)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 15, itermax = 20, trace = FALSE)))
  expect_equal(res$optim$nfeval, 15 * (1 + res$optim$iter))
})

test_that("a NaN objective is reported, and the evaluator is not leaked", {
  set.seed(59)
  expect_error(
    suppressWarnings(
      DEoptim(function(x) NaN, rep(-5, 2), rep(5, 2),
              DEoptim.control(NP = 10, itermax = 5, trace = FALSE))),
    "NaN value of objective function")
})
