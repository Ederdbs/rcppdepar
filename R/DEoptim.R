DEoptim.control <- function(VTR = -Inf, strategy = 2, bs = FALSE, NP = 50,
                            itermax = 200, CR = 0.5, F = 0.8, trace = TRUE,
                            initialpop = NULL, storepopfrom = itermax + 1,
                            storepopfreq = 1, p = 0.2, c = 0,
                            reltol = sqrt(.Machine$double.eps), steptol = itermax,
                            nthreads = 1) {
  if (itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n", immediate. = TRUE)
    itermax <- 200
  }
  ## 6, not 4: permute() draws 5 distinct indices and needs NP > 5. With NP in
  ## {4,5} the urn loop overran both of its arrays.
  if (NP < 6) {
    warning("'NP' < 6; set to default value 50\n", immediate. = TRUE)
    NP <- 50
  }
  if (F < 0 | F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n", immediate. = TRUE)
    F <- 0.8
  }
  if (CR < 0 | CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n", immediate. = TRUE)
    CR <- 0.5
  }
  if (strategy < 1 | strategy > 6) {
    warning("'strategy' not in {1,...,6}; set to default value 2\n",
            immediate. = TRUE)
    strategy <- 2
  }

  bs <- (bs > 0)

  if ( trace < 0 ) {
    warning("'trace' cannot be negative; set to 'TRUE'")
    trace <- TRUE
  }

  storepopfreq <- floor(storepopfreq)
  if (storepopfreq > itermax)
    storepopfreq <- 1
  ## storepopfreq was never checked from below; 0 reached C++ and produced a
  ## division by zero (i_iter %% i_storepopfreq), negatives were nonsense
  if (!is.finite(storepopfreq) || storepopfreq < 1) {
    warning("'storepopfreq' < 1; set to default value 1\n", immediate. = TRUE)
    storepopfreq <- 1
  }

  if (p <= 0 || p > 1) {
    warning("'p' not in (0,1]; set to default value 0.2\n", immediate. = TRUE)
    p <- 0.2
  }
  
  if (c < 0 || c > 1) {
    warning("'c' not in [0,1]; set to default value 0\n", immediate. = TRUE)
    c <- 0
  }

  if (nthreads < 1) {
    warning("'nthreads' < 1; set to default value 1\n", immediate. = TRUE)
    nthreads <- 1
  }

  list(VTR = VTR, strategy = strategy, NP = NP, itermax = itermax, CR
       = CR, F = F, bs = bs, trace = trace, initialpop = initialpop,
       storepopfrom = storepopfrom, storepopfreq = storepopfreq,
       p = p, c = c, reltol = reltol, steptol = steptol, nthreads = nthreads)
}

DEoptim <- function(fn, lower, upper, control = DEoptim.control(), ...) {
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (!is.null(names(lower)))
    nam <- names(lower)
  else if (!is.null(names(upper)) & is.null(names(lower)))
    nam <- names(upper)
  else
    nam <- paste("par", 1:length(lower), sep = "")

  env <- list2env(list(...))

  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 6) {
    warning("'NP' < 6; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  }
  else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)

  ## The parallel path is only taken for a raw-tagged compiled objective (see
  ## putFunPtrInXPtrRaw / isRawCompiledObjective in src/deoptim.cpp); anything
  ## else runs serially no matter what 'nthreads' says, so say so out loud.
  is_raw <- typeof(fn) == "externalptr" && isTRUE(attr(fn, "RcppDE::rawfun"))
  if (ctrl$nthreads > 1 && !is_raw)
    warning("'nthreads' > 1 requires a thread-safe compiled objective; ",
            "running serially on a single thread\n", immediate. = TRUE)

  ## setThreadOptions() is process-wide and used to be left changed on exit,
  ## silently altering the thread count for every other RcppParallel user in
  ## the session. Set it only when it can matter, and put it back.
  if (is_raw && ctrl$nthreads > 1) {
    old_threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA)
    on.exit(
      if (is.na(old_threads))
        Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
      else
        Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old_threads),
      add = TRUE)
    RcppParallel::setThreadOptions(numThreads = ctrl$nthreads)
  }

  outC <- DEoptim_impl(lower, upper, fn, ctrl, env)
  ##
  ## DEoptim_impl() now returns exactly the populations that were stored, so
  ## there is nothing left to slice here (the old floor() disagreed with the
  ## C++ allocation and dropped one).
  storepop <- outC$storepop
  if (length(storepop) > 0) {
    for (i in seq_along(storepop)) dimnames(storepop[[i]]) <- list(1:ctrl$NP, nam)
  } else {
    storepop <- NULL
  }

  ## optim
  bestmem <- as.numeric(outC$bestmem)
  names(bestmem) <- nam
  bestval <- as.numeric(outC$bestval)
  nfeval <- as.numeric(outC$nfeval)
  iter <- as.numeric(outC$iter)

  ## member
  names(lower) <- names(upper) <- nam
  #bestmemit <- matrix(outC$bestmemit, nrow = iter, ncol = ctrl$npar, byrow = TRUE)
  bestmemit <- outC$bestmemit[seq_len(iter), , drop = FALSE]

  ## seq_len(), not 1:iter -- iter can be 0 when VTR is already met at
  ## initialization, and 1:0 silently returns two elements
  dimnames(bestmemit) <- list(seq_len(iter), nam)
  bestvalit <- as.numeric(outC$bestvalit[seq_len(iter)])
  #pop <- matrix(outC$pop, nrow = ctrl$NP, ncol = ctrl$npar, byrow = TRUE)
  pop <- outC$pop
  storepop <- as.list(storepop)

  outR <- list(optim = list(
              bestmem = bestmem,
              bestval = bestval,
              nfeval = nfeval,
              iter = iter),
            member = list(
              lower = lower,
              upper = upper,
              bestmemit = bestmemit,
              bestvalit = bestvalit,
              pop = pop,
              storepop = storepop))

  attr(outR, "class") <- "DEoptim"
  return(outR)
}

