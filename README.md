
## RcppDEpar: Parallel Rcpp port of Differential Evolution

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.r-project.org/Licenses/GPL-2)

### About

`RcppDEpar` provides global optimization by differential evolution, with a
C++ implementation of the `DEoptim` algorithm. It is a fork of Dirk
Eddelbuettel's [RcppDE](https://github.com/eddelbuettel/rcppde) that adds a
**population-parallel evaluation path**: an entire generation's population
can be evaluated across multiple threads (via `RcppParallel`) instead of one
individual at a time, for objectives expensive enough to make that pay off.

The parallel path is opt-in and requires two things:

1. A thread-safe compiled objective, created with the internal
   `putFunPtrInXPtrRaw()` (a raw C++ function pointer, not an R closure or
   an ordinary `putFunPtrInXPtr()` pointer).
2. `control$nthreads > 1` in `DEoptim.control()`.

Everything else — plain R functions, `DEoptim.control()` options, results
format — behaves exactly like upstream `RcppDE`/`DEoptim`, and runs on the
same sequential engine when the conditions above aren't met.

See `vignette("parallel-deoptim", package = "RcppDEpar")` for a worked
example of the parallel path, including a reproducibility check (same
result regardless of thread count) and a speedup benchmark. The `demo/`
directory also has runnable scripts, including `parallel_breeding.R` for a
larger-scale parallel benchmark.

### Installation

```r
# install.packages("remotes")
remotes::install_github("Ederdbs/rcppdepar")
```

### Usage

```r
library(RcppDEpar)

Rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
fit <- DEoptim(Rosenbrock, lower = rep(-10, 2), upper = rep(10, 2))
fit$optim$bestmem
```

To use the parallel path, pass a thread-safe compiled objective (an
external pointer from `putFunPtrInXPtrRaw()`, not a plain R closure) and
set `nthreads` in `DEoptim.control()`:

```r
xptr <- RcppDEpar:::putFunPtrInXPtrRaw("breedingPenaltyRaw")  # example objective, src/exampleFunctions.cpp
fit <- DEoptim(xptr, lower = rep(0, 1000), upper = rep(1, 1000),
               control = DEoptim.control(NP = 200, nthreads = 4))
fit$optim$bestval
```

## Author

Eder David Borges da Silva, forking `RcppDE` by Dirk Eddelbuettel, which itself extends
`DEoptim` by David Ardia, Katharine Mullen, Brian Peterson and Joshua
Ulrich, based on DE-Engine by Rainer Storn.
