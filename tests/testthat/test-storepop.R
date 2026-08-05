## Regression tests for the stored-population bookkeeping.
##
## DEoptim_impl() used to size the storepop list with
##     ceil(static_cast<double>((itermax - storepopfrom) / storepopfreq))
## -- integer division, so ceil() was a no-op -- while the store branch fires
## once per generation i in [storepopfrom, itermax) with i %% storepopfreq == 0.
## Whenever storepopfreq did not divide the range evenly the list was one slot
## short and the last store wrote past the end of an R vector.
##
## The R side then re-derived the count with a *different* off-by-one
## (floor((iter - storepopfrom) / storepopfreq)) and dropped a generation.

sphere <- function(x) sum(x^2)

## Generations that store, mirroring the C++ condition, with R's 1-based
## storepopfrom converted the way DEoptim_impl() does (storepopfrom - 1).
expected_stores <- function(itermax, storepopfrom, storepopfreq) {
  gens <- seq.int(0, itermax - 1)
  sum(gens >= (storepopfrom - 1) & gens %% storepopfreq == 0)
}

run <- function(itermax, storepopfrom, storepopfreq, NP = 10) {
  set.seed(7)
  suppressWarnings(
    DEoptim(sphere, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = NP, itermax = itermax, trace = FALSE,
                            storepopfrom = storepopfrom,
                            storepopfreq = storepopfreq)))
}

test_that("storepop has exactly one entry per stored generation", {
  ## itermax=10, from=1, freq=3 is the case that used to corrupt the heap:
  ## 3 slots allocated, 4 stores (generations 0, 3, 6, 9).
  grid <- expand.grid(itermax = c(10, 11, 12), freq = c(1, 2, 3, 4, 7),
                      from = c(1, 2, 5))
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    res <- run(g$itermax, g$from, g$freq)
    expect_equal(length(res$member$storepop),
                 expected_stores(g$itermax, g$from, g$freq),
                 info = sprintf("itermax=%d from=%d freq=%d",
                                g$itermax, g$from, g$freq))
  }
})

test_that("every returned storepop entry is a populated NP x npar matrix", {
  res <- run(itermax = 10, storepopfrom = 1, storepopfreq = 3)
  expect_length(res$member$storepop, 4L)
  for (p in res$member$storepop) {
    expect_false(is.null(p))
    expect_equal(dim(p), c(10L, 2L))
    expect_true(all(is.finite(p)))
    expect_equal(colnames(p), c("par1", "par2"))
  }
})

test_that("no NULL trailing slots when the run stops before itermax", {
  ## VTR is reached long before itermax, so far fewer generations run than the
  ## list was sized for; the unfilled slots must not be handed back.
  set.seed(3)
  res <- suppressWarnings(
    DEoptim(sphere, rep(-5, 2), rep(5, 2),
            DEoptim.control(NP = 20, itermax = 500, trace = FALSE, VTR = 1e-4,
                            storepopfrom = 1, storepopfreq = 1)))
  expect_true(res$optim$iter > 0 && res$optim$iter < 500)
  expect_false(any(vapply(res$member$storepop, is.null, logical(1))))
  expect_length(res$member$storepop, res$optim$iter)
})

test_that("storepop is empty when nothing is stored (default settings)", {
  res <- run(itermax = 10, storepopfrom = 11, storepopfreq = 1)
  expect_length(res$member$storepop, 0L)
})
