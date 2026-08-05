
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2022  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#include <Rcpp/Lightest>
#include <cmath>                // pow, sin, cos, std::abs
#include <stdexcept>
#include <string>
#include "evaluate.h"

double genrose(Rcpp::NumericVector x) {       // genrose function in C++
    const double a = 1.0;
    const double b = 100.0;
    int n = x.size();
    double sum = 1.0;
    for (int i=1; i<n; i++) {
        sum += b*( ::pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - a)*(x[i] - a);
    }
    return(sum);
}

double wild(Rcpp::NumericVector x) {          // wild function in C++
    int n = x.size();
    double sum = 0.0;
    for (int i=0; i<n; i++) {
        double xsq = x[i]*x[i];
        sum += 10 * ::sin(0.3 * x[i]) * ::sin(1.3 * xsq) + 0.00001 * xsq*xsq + 0.2 * x[i] + 80;
    }
    sum /= n;
    return(sum);
}

double rastrigin(Rcpp::NumericVector x) {     // rastrigin function in C++
    int n = x.size();
    double sum = 20.0;
    for (int i=0; i<n; i++) {
        sum += x[i]+2 - 10*::cos(M_2PI*x[i]);
    }
    return(sum);
}

typedef double (*funcPtr)(Rcpp::NumericVector);

// [[Rcpp::export]]
SEXP putFunPtrInXPtr(const std::string& fstr) { 			// needed for tests/
    if (fstr == "genrose")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&genrose)));
    else if (fstr == "wild")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&wild)));
    else
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&rastrigin)));
}

// Thread-safe example objective, callable from RcppParallel worker threads
// (see EvalCompiledRaw / devol_parallel): a stand-in for a breeding-selection
// score over a candidate vector x in [0,1]^n encoding a continuous relaxation
// of "select the top-k candidates". Branchy on purpose (parent-usage cap +
// pairwise coancestry-like penalty over a local window) to represent the
// target workload: cheap per call, but not vectorizable across candidates.
double breedingPenaltyRaw(const double *x, int n) {
    const double sel_threshold = 0.5;   // candidates with x[i] > threshold are "selected"
    const int max_per_parent   = 3;     // parent-usage cap (parent = i %% 10 as a stand-in)
    const int coancestry_window = 5;    // penalize selected candidates too close together
    int parent_count[10] = {0};
    double penalty = 0.0;
    int n_selected = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] <= sel_threshold) continue;
        n_selected++;
        int parent = i % 10;
        parent_count[parent]++;
        if (parent_count[parent] > max_per_parent)
            penalty += 5.0 * (parent_count[parent] - max_per_parent);
        for (int w = 1; w <= coancestry_window && i - w >= 0; w++) {
            if (x[i - w] > sel_threshold)
                penalty += 1.0 / w;                 // closer neighbours cost more
        }
    }
    double target = n * 0.1;                        // aim to select ~10% of candidates
    penalty += std::abs(n_selected - target);
    return penalty;
}

// Raw (thread-safe) counterparts of the three example objectives above, so
// that the parallel and serial paths can be run on the *same* objective and
// compared -- see tests/testthat/test-parallel.R.
double genroseRaw(const double *x, int n) {
    const double a = 1.0, b = 100.0;
    double sum = 1.0;
    for (int i = 1; i < n; i++)
        sum += b*((x[i-1]*x[i-1] - x[i])*(x[i-1]*x[i-1] - x[i])) + (x[i] - a)*(x[i] - a);
    return sum;
}

double wildRaw(const double *x, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double xsq = x[i]*x[i];
        sum += 10 * ::sin(0.3 * x[i]) * ::sin(1.3 * xsq) + 0.00001 * xsq*xsq + 0.2 * x[i] + 80;
    }
    return sum / n;
}

double rastriginRaw(const double *x, int n) {
    double sum = 20.0;
    for (int i = 0; i < n; i++)
        sum += x[i]+2 - 10*::cos(M_2PI*x[i]);
    return sum;
}

// Always throws -- checks that a failure during the (main-thread, sequential)
// initialization loop of devol_parallel() surfaces normally.
double throwingRaw(const double *, int) {
    throw std::runtime_error("boom");
}

// Throws only after throwLateAfter calls, so that initialization (exactly NP
// evaluations) succeeds and the failure happens inside a RcppParallel worker
// thread. That is the case TinyThread's `catch(...) {}` used to swallow.
static std::atomic<long> throwLateCount(0);
static long throwLateAfter = 0;

double throwLateRaw(const double *, int) {
    if (throwLateCount.fetch_add(1, std::memory_order_relaxed) >= throwLateAfter)
        throw std::runtime_error("boom (late)");
    return 0.0;
}

// [[Rcpp::export]]
SEXP putFunPtrInXPtrRaw(const std::string& fstr, long throw_after = 0) {   // used to ignore fstr entirely
    Rcpp::DE::rawFuncPtr fun = &breedingPenaltyRaw;
    if (fstr == "genrose")        fun = &genroseRaw;
    else if (fstr == "wild")      fun = &wildRaw;
    else if (fstr == "rastrigin") fun = &rastriginRaw;
    else if (fstr == "throwing")  fun = &throwingRaw;
    else if (fstr == "throwlate") {
        fun = &throwLateRaw;
        throwLateAfter = throw_after;                // reset on each handout
        throwLateCount.store(0, std::memory_order_relaxed);
    }
    else if (fstr != "breeding" && fstr != "breedingPenaltyRaw")
        Rcpp::stop("unknown raw objective '%s'", fstr);

    Rcpp::XPtr<Rcpp::DE::rawFuncPtr> xptr(new Rcpp::DE::rawFuncPtr(fun));
    xptr.attr("RcppDE::rawfun") = true;              // tag so devol() can dispatch to the parallel path
    return xptr;
}
