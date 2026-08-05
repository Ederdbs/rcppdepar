// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Port of DEoptim (2.0.7) to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2015  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004  
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// #include <google/profiler.h>

// NB: parameter names here must stay in the same order as the definitions in
// src/devol.cpp -- they used to read (i_storepopfreq, i_storepopfrom) against a
// definition of (i_storepopfrom, i_storepopfreq). Both are int so it compiled
// and the call site happened to pass them in definition order, but the
// declaration was actively misleading.
void devol(double VTR, double f_weight, double fcross, int i_bs_flag,
           const arma::colvec & lower, const arma::colvec & upper, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax,
           arma::mat & initpopm, int i_storepopfrom, int i_storepopfreq,
           int i_specinitialpop,
           arma::mat    & ta_popP, arma::mat    & ta_oldP, arma::mat    & ta_newP, arma::colvec & t_bestP,
           arma::colvec & ta_popC, arma::colvec & ta_oldC, arma::colvec & ta_newC, double       & t_bestC,
           arma::colvec & t_bestitP, arma::colvec & t_tmpP, arma::mat & d_pop, Rcpp::List & d_storepop,
           arma::mat & d_bestmemit, arma::colvec & d_bestvalit, int & i_iterations, double i_pPct, double d_c,
           double & l_nfeval, double d_reltol, int i_steptol, int & i_storepopcnt);

// population-parallel variant, see src/devol.cpp; only usable with a
// thread-safe compiled objective (EvalCompiledRaw, tagged "RcppDE::rawfun")
void devol_parallel(double VTR, double f_weight, double f_cross, int i_bs_flag,
           const arma::colvec & fa_minbound, const arma::colvec & fa_maxbound, SEXP fcall, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax, arma::mat & initialpopm,
           int i_storepopfrom, int i_storepopfreq, int i_specinitialpop,
           arma::mat &ta_popP, arma::mat &ta_oldP, arma::mat &ta_newP, arma::colvec & t_bestP,
           arma::colvec & ta_popC, arma::colvec & ta_oldC, arma::colvec & ta_newC, double & t_bestC,
           arma::colvec & t_bestitP,
           arma::mat &d_pop, Rcpp::List &d_storepop, arma::mat & d_bestmemit, arma::colvec & d_bestvalit,
           int & i_iterations, double i_pPct, double d_c, double & l_nfeval,
           double d_reltol, int i_steptol, int i_nthreads, int & i_storepopcnt);

// Number of generations i in [i_storepopfrom, i_itermax) with i %% i_storepopfreq == 0,
// i.e. exactly how many times the store-population branch in devol()/devol_parallel()
// can fire. The old expression was
//     ceil(static_cast<double>((i_itermax - i_storepopfrom) / i_storepopfreq))
// where the division is integer division and the cast happens after it, so ceil()
// was a no-op and the result was a floor. With e.g. itermax=10, from=0, freq=3 that
// allocated 3 slots for 4 stores (i = 0,3,6,9) and Rcpp::List::operator[] does not
// bounds-check, so the fourth store corrupted the heap.
static int nStorePop(int i_itermax, int i_storepopfrom, int i_storepopfreq) {
    if (i_itermax <= 0 || i_storepopfreq < 1) return 0;
    int from = i_storepopfrom < 0 ? 0 : i_storepopfrom;
    if (from >= i_itermax) return 0;
    int first = ((from + i_storepopfreq - 1) / i_storepopfreq) * i_storepopfreq;  // first multiple >= from
    if (first >= i_itermax) return 0;
    return (i_itermax - 1 - first) / i_storepopfreq + 1;
}

// true iff fcall is an external pointer created by putFunPtrInXPtrRaw(),
// i.e. a thread-safe objective usable from RcppParallel worker threads
static bool isRawCompiledObjective(SEXP fcall) {
    if (TYPEOF(fcall) != EXTPTRSXP) return false;
    SEXP tag = Rf_getAttrib(fcall, Rf_install("RcppDE::rawfun"));
    return TYPEOF(tag) == LGLSXP && Rf_asLogical(tag) == TRUE;
}

// [[Rcpp::export]]
Rcpp::List DEoptim_impl(const arma::colvec & minbound,                  // user-defined lower bounds
                        const arma::colvec & maxbound,                  // user-defined upper bounds
                        SEXP fnS,                                       // function to be optimized, either R or C++
                        const Rcpp::List & control,                     // parameters 
                        SEXP rhoS) {                                    // optional environment
    
#if RCPP_DEV_VERSION >= RcppDevVersion(0,12,17,1)
    Rcpp::SuspendRNGSynchronizationScope rngScope;
#endif

    double VTR           = Rcpp::as<double>(control["VTR"]);            // value to reach
    int i_strategy       = Rcpp::as<int>(control["strategy"]);          // chooses DE-strategy
    int i_itermax        = Rcpp::as<int>(control["itermax"]);           // Maximum number of generations
    double l_nfeval      = 0;                                           // nb of function evals (NOT passed in)
    int i_D              = Rcpp::as<int>(control["npar"]);              // Dimension of parameter vector
    int i_NP             = Rcpp::as<int>(control["NP"]);                // Number of population members
    int i_storepopfrom   = Rcpp::as<int>(control["storepopfrom"]) - 1;  // When to start storing populations 
    int i_storepopfreq   = Rcpp::as<int>(control["storepopfreq"]);      // How often to store populations 
    int i_specinitialpop = Rcpp::as<int>(control["specinitialpop"]);    // User-defined inital population 
    double f_weight      = Rcpp::as<double>(control["F"]);              // stepsize 
    double f_cross       = Rcpp::as<double>(control["CR"]);             // crossover probability 
    int i_bs_flag        = Rcpp::as<int>(control["bs"]);                // Best of parent and child 
    int i_trace          = Rcpp::as<int>(control["trace"]);             // Print progress? 
    double i_pPct        = Rcpp::as<double>(control["p"]);              // p to define the top 100p% best solutions 
    double d_c           = Rcpp::as<double>(control["c"]);              // c as a trigger of the JADE algorithm
    double d_reltol      = Rcpp::as<double>(control["reltol"]);         // tolerance for relative convergence test, default to be sqrt(.Machine$double.eps)
    int i_steptol        = Rcpp::as<int>(control["steptol"]);           // maximum of iteration after relative convergence test is passed, default to be itermax
    int i_nthreads       = Rcpp::as<int>(control["nthreads"]);          // RcppParallel worker threads (only used with a raw compiled objective)

    // permute() draws urn_depth == 5 distinct indices and requires
    // urn_depth < NP; with NP == 4 it ran one iteration too many and wrote
    // ia_urn2[5] (array of 5) and ia_urn1[-1]. DEoptim.control() enforces this
    // too -- kept here so the C++ entry point is not overrunnable on its own.
    if (i_NP < 6) Rcpp::stop("'NP' must be at least 6");
    if (i_storepopfreq < 1) Rcpp::stop("'storepopfreq' must be at least 1");
    if (i_itermax < 1) Rcpp::stop("'itermax' must be at least 1");

    // as above, doing it in two steps is faster
    Rcpp::NumericMatrix initialpopm = Rcpp::as<Rcpp::NumericMatrix>(control["initialpop"]);
    arma::mat initpopm(initialpopm.begin(), initialpopm.rows(), initialpopm.cols(), false);

    arma::mat ta_popP(i_D, i_NP*2);                                     // Data structures for parameter vectors 
    arma::mat ta_oldP(i_D, i_NP);
    arma::mat ta_newP(i_D, i_NP);
    arma::colvec t_bestP(i_D); 

    arma::colvec ta_popC(i_NP*2);                                       // Data structures for obj. fun. values 
    arma::colvec ta_oldC(i_NP);
    arma::colvec ta_newC(i_NP);
    double t_bestC; 

    arma::colvec t_bestitP(i_D);
    arma::colvec t_tmpP(i_D); 

    int i_nstorepop = nStorePop(i_itermax, i_storepopfrom, i_storepopfreq);
    arma::mat d_pop(i_D, i_NP);
    Rcpp::List d_storepop(i_nstorepop);
    arma::mat d_bestmemit(i_D, i_itermax);
    arma::colvec d_bestvalit(i_itermax);
    int i_iter = 0;
    int i_storepopcnt = 0;                       // slots of d_storepop actually filled (the loop may stop early)

    // call actual Differential Evolution optimization given the parameters.
    // A raw-tagged compiled objective (see isRawCompiledObjective()) must
    // always go through devol_parallel()/EvalCompiledRaw, even when
    // i_nthreads == 1 (which just runs it on a single worker): its function
    // pointer has a different signature (const double*, int) than the
    // legacy EvalCompiled path expects (NumericVector), so calling it
    // through devol()/EvalCompiled would be undefined behavior.
    if (isRawCompiledObjective(fnS)) {
        devol_parallel(VTR, f_weight, f_cross, i_bs_flag, minbound, maxbound, fnS, i_trace, i_strategy, i_D, i_NP,
              i_itermax, initpopm, i_storepopfrom, i_storepopfreq, i_specinitialpop,
              ta_popP, ta_oldP, ta_newP, t_bestP, ta_popC, ta_oldC, ta_newC, t_bestC, t_bestitP,
              d_pop, d_storepop, d_bestmemit, d_bestvalit, i_iter, i_pPct, d_c, l_nfeval,
              d_reltol, i_steptol, i_nthreads, i_storepopcnt);
    } else {
        devol(VTR, f_weight, f_cross, i_bs_flag, minbound, maxbound, fnS, rhoS, i_trace, i_strategy, i_D, i_NP,
              i_itermax, initpopm, i_storepopfrom, i_storepopfreq, i_specinitialpop,
              ta_popP, ta_oldP, ta_newP, t_bestP, ta_popC, ta_oldC, ta_newC, t_bestC, t_bestitP, t_tmpP,
              d_pop, d_storepop, d_bestmemit, d_bestvalit, i_iter, i_pPct, d_c, l_nfeval,
              d_reltol, i_steptol, i_storepopcnt);
    }

    // hand back only the slots that were filled, so R does not have to
    // re-derive the count (it used to, with a different off-by-one)
    Rcpp::List storepop_out(i_storepopcnt);
    for (int i = 0; i < i_storepopcnt; i++) storepop_out[i] = d_storepop[i];

    return Rcpp::List::create(Rcpp::Named("bestmem")   = t_bestP,   // and return a named list with results to R
                              Rcpp::Named("bestval")   = t_bestC,
                              Rcpp::Named("nfeval")    = l_nfeval,
                              Rcpp::Named("iter")      = i_iter,
                              Rcpp::Named("bestmemit") = trans(d_bestmemit),
                              Rcpp::Named("bestvalit") = d_bestvalit,
                              Rcpp::Named("pop")       = trans(d_pop),
                              Rcpp::Named("storepop")  = storepop_out);

}

