demo.ParallelBreeding <- function() {

    ## Benchmark for the population-parallel path (control$nthreads > 1),
    ## using a constraint-heavy compiled objective representative of a
    ## breeding-selection problem (parent-usage cap + local coancestry
    ## penalty on a continuous relaxation of "select the top candidates"),
    ## registered in src/exampleFunctions.cpp as breedingPenaltyRaw() and
    ## exposed via RcppDEpar:::putFunPtrInXPtrRaw("breeding"). Unlike the
    ## shipped Rastrigin/Wild/Genrose benchmarks (demo/large.R), this
    ## objective is branchy rather than a trivial closed-form expression,
    ## so it better represents where population-level threading pays off.
    ##
    ## Dimension is kept well below the real 50,000-candidate use case so
    ## this demo finishes in a reasonable time; scale `n` up to see the
    ## effect grow (evaluation cost per individual is roughly linear in n).

    suppressMessages(library(RcppDEpar))

    n  <- 10000
    NP <- 2000
    itermax <- 30

    lower <- rep(0, n); upper <- rep(1, n)
    xptr <- RcppDEpar:::putFunPtrInXPtrRaw("breedingPenaltyRaw")

    runWith <- function(nthreads) {
        gc()
        set.seed(42)
        ctrl <- DEoptim.control(NP = NP, 
            itermax = itermax, trace = FALSE, 
            nthreads = nthreads)
        t <- system.time(fit <- DEoptim(xptr, 
            lower, upper, control = ctrl))[3]
        data.frame(nthreads = nthreads, seconds = as.numeric(t), bestval = fit$optim$bestval)
    }

    cat("# At", format(Sys.time()), "\n")
    cat("# n =", n, " NP =", NP, " itermax =", itermax, "\n")

    threadCounts <- c(1, 2, 8, 16, 18)
    res <- do.call(rbind, lapply(threadCounts, runWith))
    res$speedup <- res$seconds[1] / res$seconds

    print(res)
    ## bestval should match across thread counts: mutation/crossover RNG is
    ## seeded per (generation, individual index), not per thread, so results
    ## are reproducible regardless of how many threads run the population loop
    stopifnot(all(abs(res$bestval - res$bestval[1]) < 1e-8))
    cat("# bestval identical across thread counts: OK\n")
    cat("# Done", format(Sys.time()), "\n")
}


demo.ParallelBreeding()
