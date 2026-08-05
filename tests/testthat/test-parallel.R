## The population-parallel path (devol_parallel), reached only for a
## raw-tagged compiled objective built by putFunPtrInXPtrRaw().

raw_fn <- function(name) RcppDEpar:::putFunPtrInXPtrRaw(name)

test_that("putFunPtrInXPtrRaw honours its argument", {
  ## It used to ignore fstr and always return breedingPenaltyRaw.
  ## genrose has a known minimum of 1 at the all-ones vector; the breeding
  ## penalty on the same box is nowhere near that.
  set.seed(21)
  res <- suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(0.5, 3), rep(1.5, 3),
            DEoptim.control(NP = 40, itermax = 200, trace = FALSE)))
  expect_equal(res$optim$bestval, 1, tolerance = 1e-4)
  expect_equal(unname(res$optim$bestmem), rep(1, 3), tolerance = 1e-2)

  expect_error(raw_fn("no-such-function"), "unknown raw objective")
})

test_that("a raw objective is tagged and dispatches to the parallel path", {
  p <- raw_fn("genrose")
  expect_identical(typeof(p), "externalptr")
  expect_true(isTRUE(attr(p, "RcppDE::rawfun")))
})

test_that("the parallel path does not overwrite a population member either", {
  ## Same unsafe_col aliasing as devol(); here t_bestP aliased a column of
  ## ta_oldP, which is then swapped with ta_newP every generation, so the
  ## writes landed in whichever buffer happened to own that memory.
  set.seed(41)
  pop0 <- matrix(runif(20 * 3, -1, 1), nrow = 20, ncol = 3)
  res <- suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(-5, 3), rep(5, 3),
            DEoptim.control(NP = 20, itermax = 1, trace = FALSE, nthreads = 2,
                            initialpop = pop0, storepopfrom = 1)))
  expect_equal(unname(res$member$storepop[[1]]), pop0, tolerance = 1e-12)
})

test_that("results are identical regardless of thread count", {
  ## Each individual seeds its own SplitMix64 stream from (generation seed,
  ## index), not from a thread id, so the answer must not depend on how
  ## RcppParallel chunks the population.
  runs <- lapply(c(1, 2, 4), function(nth) {
    set.seed(99)
    suppressWarnings(
      DEoptim(raw_fn("genrose"), rep(-2, 4), rep(2, 4),
              DEoptim.control(NP = 32, itermax = 60, trace = FALSE,
                              nthreads = nth)))
  })
  expect_equal(runs[[1]]$optim$bestval, runs[[2]]$optim$bestval)
  expect_equal(runs[[1]]$optim$bestval, runs[[3]]$optim$bestval)
  expect_equal(runs[[1]]$optim$bestmem, runs[[2]]$optim$bestmem)
  expect_equal(runs[[1]]$member$pop, runs[[3]]$member$pop)
})

test_that("the parallel path is reproducible under set.seed", {
  once <- function() {
    set.seed(1234)
    suppressWarnings(
      DEoptim(raw_fn("genrose"), rep(-2, 3), rep(2, 3),
              DEoptim.control(NP = 20, itermax = 30, trace = FALSE,
                              nthreads = 2)))
  }
  expect_equal(once()$optim$bestval, once()$optim$bestval)
  expect_equal(once()$member$bestmemit, once()$member$bestmemit)
})

test_that("an objective that throws inside a worker thread becomes an R error", {
  ## RcppParallel's TinyThread backend wraps each task in `catch(...) {}`, so
  ## before the fix the affected slice of ta_newP/ta_newC was left holding
  ## stale (generation 1: uninitialized) values and the run reported a best.
  ## throw_after = NP lets the sequential initialization succeed, so the first
  ## failure happens in a worker.
  expect_error(
    suppressWarnings(
      DEoptim(RcppDEpar:::putFunPtrInXPtrRaw("throwlate", throw_after = 20),
              rep(-2, 3), rep(2, 3),
              DEoptim.control(NP = 20, itermax = 10, trace = FALSE,
                              nthreads = 2))),
    "objective function failed for population member")
})

test_that("an objective that throws during initialization also errors", {
  ## The init loop runs on the main thread, so the C++ exception propagates
  ## normally -- what matters is that it is not swallowed.
  expect_error(
    suppressWarnings(
      DEoptim(raw_fn("throwing"), rep(-2, 3), rep(2, 3),
              DEoptim.control(NP = 20, itermax = 10, trace = FALSE,
                              nthreads = 2))),
    "boom")
})

test_that("nthreads does not leak into the rest of the session", {
  before <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA)
  set.seed(2)
  suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(-2, 3), rep(2, 3),
            DEoptim.control(NP = 20, itermax = 10, trace = FALSE, nthreads = 4)))
  after <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA)
  expect_identical(before, after)
})

test_that("nfeval is exact and does not depend on the platform's 'long'", {
  ## NP evaluations to initialize, then NP per generation.
  set.seed(8)
  ctrl <- DEoptim.control(NP = 20, itermax = 15, trace = FALSE, nthreads = 2)
  res <- suppressWarnings(DEoptim(raw_fn("genrose"), rep(-2, 3), rep(2, 3), ctrl))
  expect_type(res$optim$nfeval, "double")
  expect_equal(res$optim$nfeval, 20 * (1 + res$optim$iter))
})

test_that("serial and parallel paths agree on the optimum they find", {
  ## Different RNG streams, so not bit-identical -- but they must land on the
  ## same minimum of the same function.
  genrose_r <- function(x) 1 + sum(100 * (x[-length(x)]^2 - x[-1])^2 + (x[-1] - 1)^2)

  set.seed(17)
  ser <- suppressWarnings(
    DEoptim(genrose_r, rep(-2, 4), rep(2, 4),
            DEoptim.control(NP = 50, itermax = 400, trace = FALSE)))
  set.seed(17)
  par <- suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(-2, 4), rep(2, 4),
            DEoptim.control(NP = 50, itermax = 400, trace = FALSE, nthreads = 2)))

  expect_equal(ser$optim$bestval, 1, tolerance = 1e-3)
  expect_equal(par$optim$bestval, 1, tolerance = 1e-3)
  ## genrose only ever sees x[1] squared, so x[1] = +-1 are both optimal;
  ## compare the coordinates that are actually determined.
  expect_equal(abs(unname(par$optim$bestmem)), rep(1, 4), tolerance = 1e-2)
  expect_equal(unname(par$optim$bestmem)[-1], unname(ser$optim$bestmem)[-1],
               tolerance = 1e-2)
})

test_that("bs = TRUE works on the parallel path", {
  ## The best-of-parent-and-child branch now reuses the caller's 2*NP scratch
  ## instead of allocating join_rows/join_cols temporaries per generation.
  set.seed(31)
  res <- suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(-2, 3), rep(2, 3),
            DEoptim.control(NP = 40, itermax = 500, trace = FALSE,
                            bs = TRUE, nthreads = 2)))
  expect_equal(res$optim$bestval, 1, tolerance = 1e-3)
  expect_true(all(is.finite(res$member$pop)))
  ## bestvalit must still be monotone with best-of-both selection
  expect_true(all(diff(res$member$bestvalit) <= 0))
})

test_that("storepop bookkeeping matches on the parallel path too", {
  set.seed(4)
  res <- suppressWarnings(
    DEoptim(raw_fn("genrose"), rep(-2, 2), rep(2, 2),
            DEoptim.control(NP = 10, itermax = 10, trace = FALSE, nthreads = 2,
                            storepopfrom = 1, storepopfreq = 3)))
  expect_length(res$member$storepop, 4L)   # generations 0, 3, 6, 9
  expect_true(all(vapply(res$member$storepop,
                         function(p) all(is.finite(p)), logical(1))))
})
