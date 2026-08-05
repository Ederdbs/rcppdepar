
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2022  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// Use RcppParallel's portable TinyThread backend instead of TBB: it's
// header-only (no extra libs/rpaths to link against, unlike TBB which
// requires linking RcppParallel's own compiled glue), and the population
// loop here is not fine-grained enough to need TBB's work-stealing scheduler.
#define RCPP_PARALLEL_USE_TBB 0
#include <RcppParallel.h>
#include <random>
#include <cstdint>
#include <cstring>              // memmove
#include <cmath>                // round, pow, std::abs
#include <memory>               // unique_ptr
#include <string>
#include "evaluate.h"           // simple function evaluation framework
// #include <google/profiler.h>

void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urntmp[]);

void devol(double VTR, double f_weight, double f_cross, int i_bs_flag,
           const arma::colvec & fa_minbound, const arma::colvec & fa_maxbound, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax, arma::mat & initialpopm,
           int i_storepopfrom, int i_storepopfreq, int i_specinitialpop,
           arma::mat &ta_popP, arma::mat &ta_oldP, arma::mat &ta_newP, arma::colvec & t_bestP,
           arma::colvec & ta_popC, arma::colvec & ta_oldC, arma::colvec & ta_newC, double & t_bestC,
           arma::colvec & t_bestitP, arma::colvec & t_tmpP,
           arma::mat &d_pop, Rcpp::List &d_storepop, arma::mat & d_bestmemit, arma::colvec & d_bestvalit,
           int & i_iterations, double i_pPct, double d_c, double & l_nfeval,
           double d_reltol, int i_steptol, int & i_storepopcnt) {

    //ProfilerStart("/tmp/RcppDE.prof");
    // unique_ptr: EvalStandard::eval() calls Rcpp::stop() on a NaN objective,
    // which unwinds straight out of devol() -- a raw new here leaked
    std::unique_ptr<Rcpp::DE::EvalBase> ev;
    if (TYPEOF(fcall) == EXTPTRSXP) {                // non-standard mode: we are being passed an external pointer
        ev.reset(new Rcpp::DE::EvalCompiled(fcall, rho)); // so assign a pointer using external pointer in fcall SEXP
    } else {                                         // standard mode: env_ is an env, fcall_ is a function
        ev.reset(new Rcpp::DE::EvalStandard(fcall, rho)); // so assign R function and environment
    }
    const int urn_depth = 5;                    // 4 + one index to avoid
    Rcpp::NumericVector par(i_D);               // initialize parameter vector to pass to evaluate function
    std::vector<int32_t> ia_urn2(urn_depth);    // fixed-size vector for urn draws
    std::vector<int32_t> ia_urntmp(i_NP);       // so that we don't need to re-allocated each time in permute
    int p_NP = round(i_pPct * i_NP);            // choose at least two best solutions
    p_NP = p_NP < 2 ? 2 : p_NP;
    std::vector<int32_t> sortIndex;             // sorted-index scratch for strategy 6 only
    if (i_strategy == 6) {
        sortIndex.resize(i_NP);
        for (int i = 0; i < i_NP; i++)
            sortIndex[i] = i;
    }

    d_bestmemit.zeros();                        // initialize best members
    d_bestvalit.zeros();                        // initialize best values
    d_pop.zeros();                              // initialize best population

    arma::mat initialpop;                       // only materialized when a starting population was supplied
    if (i_specinitialpop > 0) {                 // if initial population provided, initialize with values
        initialpop = trans(initialpopm);        // transpose as we prefer columns for population members here
    }

    for (int i = 0; i < i_NP; i++) {            // ------Initialization-----------------------------
        if (i_specinitialpop <= 0) {            // random initial member
            for (int j = 0; j < i_D; j++) {
                ta_popP.at(j,i) = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
            }
        } else {                                // or user-specified initial member
            ta_popP.col(i) = initialpop.col(i);
        }
        memmove(REAL(par), ta_popP.colptr(i), Rf_nrows(par) * sizeof(double));
        ta_popC[i] = ev->eval(par);
        if (i == 0 || ta_popC[i] <= t_bestC) {
            t_bestC = ta_popC[i];
            // .col(), not .unsafe_col(): assigning an unsafe_col to a colvec
            // makes the colvec an *alias* of that population column (arma
            // adopts the external memory rather than copying into it). The
            // next best-update then wrote the new best member straight into
            // the previously-aliased column, silently destroying a population
            // member -- always member 0 at minimum, since i == 0 aliases first.
            t_bestP = ta_popP.col(i);
        }
    }

    ta_oldP = ta_popP.cols(0, i_NP-1);          // ---assign pointers to current ("old") population---
    ta_oldC = ta_popC.rows(0, i_NP-1);

    int i_iter = 0;                             // ------Iteration loop--------------------------------------------
    int popcnt = 0;
    int i_iter_tol = 0;

    //Trigger JADE algorithm when d_c <> 0 (randomize cross-over and weight coefficient)
    //Using user supplied cross-over and weight coefficient as the mean of randomized coefficients

    double meanCR = f_cross, meanF = f_weight;

    double rand_cross  = f_cross;
    double rand_weight = f_weight;

    while ((i_iter < i_itermax) && (t_bestC > VTR) && (i_iter_tol <= i_steptol)) {    // main loop ====================================
        Rcpp::checkUserInterrupt();             // throws (does not longjmp), so C++ destructors still run

        if (i_iter % i_storepopfreq == 0 && i_iter >= i_storepopfrom) {         // store intermediate populations
            d_storepop[popcnt++] = Rcpp::wrap( trans(ta_oldP) );
        } // end store pop

        d_bestmemit.col(i_iter) = t_bestP;      // store the best member
        d_bestvalit[i_iter] = t_bestC;          // store the best value
        t_bestitP = t_bestP;
        i_iter++;                               // increase iteration counter

        double f_dither = f_weight + ::unif_rand() * (1.0 - f_weight);  // ----computer dithering factor --------------

        if (i_strategy == 6) {                  // ---DE/current-to-p-best/1 ------------------------------------------
            arma::uvec ord = arma::sort_index(ta_oldC);                 // ascending order, no copy of ta_oldC needed
            std::copy(ord.begin(), ord.end(), sortIndex.begin());
        }

        // JADE success statistics: reset every generation. These used to be
        // declared outside the while loop and never zeroed, so from the second
        // generation on the adaptation was driven by all-time accumulations
        // rather than this generation's successful trials.
        double sumCR = 0.0, sumF = 0.0, sumF2 = 0.0;
        int nGood = 0;

        for (int i = 0; i < i_NP; i++) {        // ----start of loop through ensemble------------------------

            t_tmpP = ta_oldP.col(i);            // t_tmpP is the vector to mutate and eventually select

            permute(ia_urn2.data(), urn_depth, i_NP, i, ia_urntmp.data()); // Pick 4 random and distinct

            //Trigger JADE algorithm when d_c <> 0 (randomize cross-over and weight coefficient)
            if(d_c>0){
              rand_cross = R::rnorm(meanCR, 0.1);//changed from f_cross to meanCR
              rand_cross = rand_cross > 1.0 ? 1 : rand_cross;
              rand_cross = rand_cross < 0.0 ? 0 : rand_cross;
              do{
                rand_weight = R::rcauchy(meanF, 0.1);//changed from f_weight to meanF
                rand_weight = rand_weight > 1.0 ? 1.0 : rand_weight;
              } while (rand_weight <= 0.0);
            }

            int k = 0;                          // loop counter used in all strategies below

            // ===Choice of strategy=======================================================
            switch (i_strategy) {

            case 1: {                           // ---classical strategy DE/rand/1/bin---------------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target
                    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + rand_weight *
                        (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            case 2: {                           // ---DE/local-to-best/1/bin-------------------------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target
                    t_tmpP[j] = t_tmpP[j] + rand_weight * (t_bestitP[j] - t_tmpP[j]) + rand_weight *
                        (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            case 3: {                           // ---DE/best/1/bin with jitter---------------------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target
                    double f_jitter = 0.0001 * ::unif_rand() + rand_weight;
                    t_tmpP[j] = t_bestitP[j] + f_jitter * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            case 4: {                           // ---DE/rand/1/bin with per-vector-dither----------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target *
                    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + (rand_weight + ::unif_rand()*(1.0 - rand_weight))
                        * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            case 5: {                           // ---DE/rand/1/bin with per-generation-dither------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target
                    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + f_dither
                        * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            case 6: {                           // ---DE/current-to-p-best/1 (JADE)-----------------------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                do {                            // add fluctuation to random target
                    int i_pbest = sortIndex[static_cast<int>(::unif_rand() * p_NP)]; // select from [0, 1, 2, ..., (pNP-1)]
                    t_tmpP[j] = ta_oldP.at(j,i) + rand_weight * (ta_oldP.at(j,i_pbest) - ta_oldP.at(j,i)) +
                        rand_weight * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
                    j = (j + 1) % i_D;
                } while ((::unif_rand() < rand_cross) && (++k < i_D));
                break;
            }
            default: {                          // ---variation to DE/rand/1/bin: either-or-algorithm------------------
                int j = static_cast<int>(::unif_rand() * i_D);  // random parameter
                if (::unif_rand() < 0.5) {      // differential mutation, Pmu = 0.5
                    do {                        // add fluctuation to random target */
                        t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + rand_weight * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                        j = (j + 1) % i_D;
                    } while ((::unif_rand() < rand_cross) && (++k < i_D));

                } else {                        // recombination with K = 0.5*(F+1) -. F-K-Rule
                    do {                        // add fluctuation to random target */
                        t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + 0.5 * (rand_weight + 1.0) *
                            (ta_oldP.at(j,ia_urn2[2]) + ta_oldP.at(j,ia_urn2[3]) - 2 * ta_oldP.at(j,ia_urn2[1]));
                        j = (j + 1) % i_D;
                    } while ((::unif_rand() < rand_cross) && (++k < i_D));
                }
                break;
            }
            } // end switch (i_strategy) ...

            for (int j = 0; j < i_D; j++) {     // boundary constraints, bounce-back method was not enforcing bounds
                if (t_tmpP[j] < fa_minbound[j]) {
                    t_tmpP[j] = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
                }
                if (t_tmpP[j] > fa_maxbound[j]) {
                    t_tmpP[j] = fa_maxbound[j] - ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
                }
            }

            // ------Trial mutation now in t_tmpP-----------------
            memmove(REAL(par), t_tmpP.memptr(), Rf_nrows(par) * sizeof(double));
            double t_tmpC = ev->eval(par);                              // Evaluate mutant in t_tmpP
            if (t_tmpC <= ta_oldC[i] || i_bs_flag) {                    // i_bs_flag means will choose best NP later
                ta_newP.col(i) = t_tmpP;                                // replace target with mutant
                ta_newC[i] = t_tmpC;
                if (t_tmpC <= t_bestC) {
                    t_bestP = t_tmpP;
                    t_bestC = t_tmpC;
                }
                if (d_c > 0) {                  // accumulate this generation's successful CR / F
                    sumCR += rand_cross;
                    sumF  += rand_weight;
                    sumF2 += rand_weight * rand_weight;
                    nGood++;
                }
            } else {
                ta_newP.col(i) = ta_oldP.col(i);
                ta_newC[i] = ta_oldC[i];
            }
        } // End mutation loop through pop., ie the "for (i = 0; i < i_NP; i++)"

        // guarded on nGood: with no successful trial this was 0.0/0.0 = NaN,
        // which poisoned meanF and then every subsequent trial vector
        if (d_c > 0 && nGood > 0) {             /* calculate new meanCR and meanF */
            meanCR = (1-d_c)*meanCR + d_c*(sumCR/nGood);        // arithmetic mean of successful CR
            meanF  = (1-d_c)*meanF  + d_c*(sumF2/sumF);         // Lehmer mean of successful F
        }

        if (i_bs_flag) {        // examine old and new pop. and take the best NP members into next generation

            ta_popP.cols(0, i_NP-1) = ta_oldP;
            ta_popC.rows(0, i_NP-1) = ta_oldC;

            ta_popP.cols(i_NP, 2*i_NP-1) = ta_newP;
            ta_popC.rows(i_NP, 2*i_NP-1) = ta_newC;

            arma::uvec ord = arma::sort_index(ta_popC);  // O(n log n) instead of the previous O(n^2) comb sort
            for (int i = 0; i < i_NP; i++) {              // keep only the best NP members, in ascending order
                ta_newP.col(i) = ta_popP.col(ord[i]);
                ta_newC[i] = ta_popC[ord[i]];
            }
        } // i_bs_flag

        ta_oldP.swap(ta_newP);                  // have selected NP mutants move on to next generation
        ta_oldC.swap(ta_newC);                  // swap, not copy: both are owned by the caller and unaliased

        // print out temporary results
        if ( (i_trace > 0)  &&  ((i_iter % i_trace) == 0) ) {
            Rprintf("Iteration: %d/%d bestvalit: %f\n", i_iter, i_itermax, t_bestC);
        }

        //check relative convergence
        if (std::abs(d_bestvalit[i_iter - 1] - t_bestC) < (d_reltol*(std::abs(d_bestvalit[i_iter - 1]) + d_reltol))) {
            i_iter_tol++;
        } else {
            i_iter_tol = 0;
        }

    } // end loop through generations

    d_pop = ta_oldP;
    i_iterations = i_iter;
    i_storepopcnt = popcnt;                     // how many slots of d_storepop were actually filled
    l_nfeval = static_cast<double>(ev->getNbEvals());
    // ProfilerStop();
}

// ponytail: devol_parallel() duplicates the mutation/crossover switch from
// devol() above instead of sharing it. The two loops need genuinely
// different RNG sources (R's ::unif_rand()/R::rnorm/R::rcauchy are not
// thread-safe; the parallel path needs a per-index splitmix64 instead),
// and templating the switch on an RNG policy risked subtly changing the
// well-tested serial path for a very small code-size win. If a third RNG
// backend is ever needed, factor the switch into a template then.
namespace {

    // splitmix64 (Vigna): single 64-bit word of state, O(1) to seed, unlike
    // std::mt19937_64's ~312-word state (re-initialized from scratch for
    // every individual below). Satisfies UniformRandomBitGenerator, so it
    // drops straight into std::uniform_real_distribution/normal/cauchy.
    struct SplitMix64 {
        using result_type = uint64_t;
        uint64_t state;
        explicit SplitMix64(uint64_t seed) : state(seed) {}
        static constexpr result_type min() { return 0; }
        static constexpr result_type max() { return UINT64_MAX; }
        result_type operator()() {
            uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            return z ^ (z >> 31);
        }
    };

    // Deterministic per-individual RNG: seeded from (generation base seed, index),
    // not from thread id, so results don't depend on how RcppParallel schedules
    // work across threads/thread counts.
    struct ThreadRng {
        SplitMix64 gen;
        explicit ThreadRng(uint64_t seed) : gen(seed) {}
        double unif()                 { return std::uniform_real_distribution<double>(0.0, 1.0)(gen); }
        double norm(double m, double s)   { return std::normal_distribution<double>(m, s)(gen); }
        double cauchy(double m, double s) { return std::cauchy_distribution<double>(m, s)(gen); }
    };

    // Thread-safe re-implementation of permute() (src/permute.cpp), drawing
    // from a local ThreadRng instead of R's global ::unif_rand().
    void permuteTS(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urn1[], ThreadRng &rng) {
        int k = i_NP;
        int i_urn1 = 0;
        int i_urn2 = 0;
        for (int i = 0; i < i_NP; i++)
            ia_urn1[i] = i;
        i_urn1 = i_avoid;
        while (k > i_NP - i_urn2_depth) {
            ia_urn2[i_urn2] = ia_urn1[i_urn1];
            ia_urn1[i_urn1] = ia_urn1[k-1];
            k = k - 1;
            i_urn2 = i_urn2 + 1;
            i_urn1 = static_cast<int>(rng.unif() * k);
        }
    }

    struct DEWorker : public RcppParallel::Worker {
        // ---- read-only shared inputs ----
        const arma::mat &ta_oldP;
        const arma::colvec &ta_oldC;
        const arma::colvec &t_bestitP;
        const std::vector<int32_t> &sortIndex;
        const arma::colvec &fa_minbound, &fa_maxbound;
        int i_D, i_NP, i_strategy, i_bs_flag, p_NP;
        double f_weight, f_cross, f_dither, d_c;
        uint64_t baseSeed;
        Rcpp::DE::EvalCompiledRaw *ev;

        // ---- outputs: index i is written only by the task handling i ----
        arma::mat &ta_newP;
        arma::colvec &ta_newC;
        arma::colvec &out_rand_cross;
        arma::colvec &out_rand_weight;
        std::vector<unsigned char> &accepted;

        // RcppParallel's TinyThread backend wraps each task in `catch(...) {}`
        // and swallows the exception, which would leave this task's slice of
        // ta_newP/ta_newC untouched (stale or, in generation 1, uninitialized)
        // while the optimization carried on reporting a bogus best. Record the
        // first failing index so the main thread can turn it into an R error.
        std::atomic<int> errIndex;

        DEWorker(const arma::mat &ta_oldP_, const arma::colvec &ta_oldC_, const arma::colvec &t_bestitP_,
                 const std::vector<int32_t> &sortIndex_, const arma::colvec &fa_minbound_, const arma::colvec &fa_maxbound_,
                 int i_D_, int i_NP_, int i_strategy_, int i_bs_flag_, int p_NP_,
                 double f_weight_, double f_cross_, double f_dither_, double d_c_, uint64_t baseSeed_,
                 Rcpp::DE::EvalCompiledRaw *ev_,
                 arma::mat &ta_newP_, arma::colvec &ta_newC_,
                 arma::colvec &out_rand_cross_, arma::colvec &out_rand_weight_, std::vector<unsigned char> &accepted_)
            : ta_oldP(ta_oldP_), ta_oldC(ta_oldC_), t_bestitP(t_bestitP_), sortIndex(sortIndex_),
              fa_minbound(fa_minbound_), fa_maxbound(fa_maxbound_),
              i_D(i_D_), i_NP(i_NP_), i_strategy(i_strategy_), i_bs_flag(i_bs_flag_), p_NP(p_NP_),
              f_weight(f_weight_), f_cross(f_cross_), f_dither(f_dither_), d_c(d_c_), baseSeed(baseSeed_), ev(ev_),
              ta_newP(ta_newP_), ta_newC(ta_newC_),
              out_rand_cross(out_rand_cross_), out_rand_weight(out_rand_weight_), accepted(accepted_),
              errIndex(-1) {}

        void operator()(std::size_t begin, std::size_t end) {
            std::vector<int32_t> ia_urn2(5);
            std::vector<int32_t> ia_urntmp(i_NP);
            arma::colvec t_tmpP(i_D);

            for (std::size_t idx = begin; idx < end; idx++) {
                int i = static_cast<int>(idx);
                try {
                ThreadRng rng(baseSeed ^ (0x9E3779B97F4A7C15ULL * static_cast<uint64_t>(i + 1)));

                t_tmpP = ta_oldP.col(i);
                permuteTS(ia_urn2.data(), 5, i_NP, i, ia_urntmp.data(), rng);

                double rand_cross = f_cross, rand_weight = f_weight;
                if (d_c > 0) {
                    rand_cross = rng.norm(f_cross, 0.1);
                    rand_cross = rand_cross > 1.0 ? 1 : rand_cross;
                    rand_cross = rand_cross < 0.0 ? 0 : rand_cross;
                    do {
                        rand_weight = rng.cauchy(f_weight, 0.1);
                        rand_weight = rand_weight > 1.0 ? 1.0 : rand_weight;
                    } while (rand_weight <= 0.0);
                }

                int k = 0;
                int j = static_cast<int>(rng.unif() * i_D);

                switch (i_strategy) {
                case 1: {
                    do {
                        t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + rand_weight *
                            (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                case 2: {
                    do {
                        t_tmpP[j] = t_tmpP[j] + rand_weight * (t_bestitP[j] - t_tmpP[j]) + rand_weight *
                            (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                case 3: {
                    do {
                        double f_jitter = 0.0001 * rng.unif() + rand_weight;
                        t_tmpP[j] = t_bestitP[j] + f_jitter * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                case 4: {
                    do {
                        t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + (rand_weight + rng.unif()*(1.0 - rand_weight))
                            * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                case 5: {
                    do {
                        t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + f_dither
                            * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                case 6: {
                    do {
                        int i_pbest = sortIndex[static_cast<int>(rng.unif() * p_NP)];
                        t_tmpP[j] = ta_oldP.at(j,i) + rand_weight * (ta_oldP.at(j,i_pbest) - ta_oldP.at(j,i)) +
                            rand_weight * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
                        j = (j + 1) % i_D;
                    } while ((rng.unif() < rand_cross) && (++k < i_D));
                    break;
                }
                default: {
                    if (rng.unif() < 0.5) {
                        do {
                            t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + rand_weight * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
                            j = (j + 1) % i_D;
                        } while ((rng.unif() < rand_cross) && (++k < i_D));
                    } else {
                        do {
                            t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + 0.5 * (rand_weight + 1.0) *
                                (ta_oldP.at(j,ia_urn2[2]) + ta_oldP.at(j,ia_urn2[3]) - 2 * ta_oldP.at(j,ia_urn2[1]));
                            j = (j + 1) % i_D;
                        } while ((rng.unif() < rand_cross) && (++k < i_D));
                    }
                    break;
                }
                } // end switch

                for (int jj = 0; jj < i_D; jj++) {
                    if (t_tmpP[jj] < fa_minbound[jj])
                        t_tmpP[jj] = fa_minbound[jj] + rng.unif() * (fa_maxbound[jj] - fa_minbound[jj]);
                    if (t_tmpP[jj] > fa_maxbound[jj])
                        t_tmpP[jj] = fa_maxbound[jj] - rng.unif() * (fa_maxbound[jj] - fa_minbound[jj]);
                }

                double t_tmpC = ev->evalRaw(t_tmpP.memptr(), i_D);   // thread-safe: no R API calls
                if (t_tmpC <= ta_oldC[i] || i_bs_flag) {
                    ta_newP.col(i) = t_tmpP;
                    ta_newC[i] = t_tmpC;
                    accepted[i] = 1;
                    out_rand_cross[i] = rand_cross;
                    out_rand_weight[i] = rand_weight;
                } else {
                    ta_newP.col(i) = ta_oldP.col(i);
                    ta_newC[i] = ta_oldC[i];
                    accepted[i] = 0;
                }
                } catch (...) {                 // see errIndex above; record the first failure and bail out
                    int expected = -1;
                    errIndex.compare_exchange_strong(expected, i);
                    return;
                }
            }
        }
    };

} // anonymous namespace

// Population-parallel variant of devol(), dispatched from DEoptim_impl only
// when the objective is a thread-safe compiled function (EvalCompiledRaw,
// tagged "RcppDE::rawfun" by putFunPtrInXPtrRaw()) and control$nthreads > 1.
// The generation loop stays sequential (each generation depends on the last);
// only the per-individual mutation/evaluation loop runs across threads.
//
// Best-tracking and JADE accumulation are moved to a post-loop reduction
// instead of updating shared state inside the loop. This is safe because
// t_bestitP (the only best-so-far value read during mutation) is snapshotted
// once at the top of each generation (mirroring devol()) -- intra-generation
// updates to t_bestP never feed back into the same generation's mutations in
// either implementation, so deferring the update changes nothing observable.
void devol_parallel(double VTR, double f_weight, double f_cross, int i_bs_flag,
           const arma::colvec & fa_minbound, const arma::colvec & fa_maxbound, SEXP fcall, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax, arma::mat & initialpopm,
           int i_storepopfrom, int i_storepopfreq, int i_specinitialpop,
           arma::mat &ta_popP, arma::mat &ta_oldP, arma::mat &ta_newP, arma::colvec & t_bestP,
           arma::colvec & ta_popC, arma::colvec & ta_oldC, arma::colvec & ta_newC, double & t_bestC,
           arma::colvec & t_bestitP,
           arma::mat &d_pop, Rcpp::List &d_storepop, arma::mat & d_bestmemit, arma::colvec & d_bestvalit,
           int & i_iterations, double i_pPct, double d_c, double & l_nfeval,
           double d_reltol, int i_steptol, int i_nthreads, int & i_storepopcnt) {

    // thread count is set at the R level via RcppParallel::setThreadOptions()
    // (see R/DEoptim.R) before this function is entered; i_nthreads is kept
    // as a parameter for clarity/future use but not read here
    (void) i_nthreads;

    Rcpp::DE::EvalCompiledRaw ev(fcall);        // caller guarantees fcall is raw-tagged (thread-safe)

    int p_NP = round(i_pPct * i_NP);
    p_NP = p_NP < 2 ? 2 : p_NP;
    std::vector<int32_t> sortIndex;             // sorted-index scratch for strategy 6 only
    if (i_strategy == 6) {
        sortIndex.resize(i_NP);
        for (int i = 0; i < i_NP; i++)
            sortIndex[i] = i;
    }

    d_bestmemit.zeros();
    d_bestvalit.zeros();
    d_pop.zeros();

    arma::mat initialpop;                       // only materialized when a starting population was supplied
    if (i_specinitialpop > 0) {
        initialpop = trans(initialpopm);
    }

    arma::colvec initC(i_NP);
    for (int i = 0; i < i_NP; i++) {                          // ------Initialization----- (sequential, main thread)
        if (i_specinitialpop <= 0) {
            for (int j = 0; j < i_D; j++) {
                ta_oldP.at(j,i) = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
            }
        } else {
            ta_oldP.col(i) = initialpop.col(i);
        }
        initC[i] = ev.evalRaw(ta_oldP.colptr(i), i_D);
        if (i == 0 || initC[i] <= t_bestC) {
            t_bestC = initC[i];
            t_bestP = ta_oldP.col(i);       // .col(), not .unsafe_col() -- see devol() above
        }
    }
    ta_oldC = initC;

    int i_iter = 0;
    int popcnt = 0;
    int i_iter_tol = 0;

    double meanCR = f_cross, meanF = f_weight;

    arma::colvec rand_cross_used(i_NP), rand_weight_used(i_NP);
    std::vector<unsigned char> accepted(i_NP);

    while ((i_iter < i_itermax) && (t_bestC > VTR) && (i_iter_tol <= i_steptol)) {
        Rcpp::checkUserInterrupt();             // main thread only, between generations

        if (i_iter % i_storepopfreq == 0 && i_iter >= i_storepopfrom) {
            d_storepop[popcnt++] = Rcpp::wrap( trans(ta_oldP) );
        }

        d_bestmemit.col(i_iter) = t_bestP;
        d_bestvalit[i_iter] = t_bestC;
        t_bestitP = t_bestP;
        i_iter++;

        double f_dither = f_weight + ::unif_rand() * (1.0 - f_weight);

        if (i_strategy == 6) {
            arma::uvec ord = arma::sort_index(ta_oldC);
            std::copy(ord.begin(), ord.end(), sortIndex.begin());
        }

        // one base seed drawn from R's RNG (main thread, safe) per generation;
        // each individual then derives its own deterministic stream from it.
        // Scale by 2^53, not 2^64: 18446744073709551615.0 is not representable
        // as a double and rounds to exactly 2^64, so a unif_rand() close enough
        // to 1.0 made the float->uint64_t conversion undefined (and the low 12
        // bits of the seed were always zero). SplitMix64 diffuses the 53 bits.
        uint64_t baseSeed = static_cast<uint64_t>(::unif_rand() * 9007199254740992.0);

        // when d_c>0 (JADE), mutation draws CR/F around the running means
        // meanCR/meanF rather than the fixed f_cross/f_weight, matching
        // devol()'s use of meanCR/meanF at src/devol.cpp:114,118
        double jadeCR = (d_c > 0) ? meanCR : f_cross;
        double jadeF  = (d_c > 0) ? meanF  : f_weight;

        DEWorker worker(ta_oldP, ta_oldC, t_bestitP, sortIndex, fa_minbound, fa_maxbound,
                         i_D, i_NP, i_strategy, i_bs_flag, p_NP,
                         jadeF, jadeCR, f_dither, d_c, baseSeed, &ev,
                         ta_newP, ta_newC, rand_cross_used, rand_weight_used, accepted);

        RcppParallel::parallelFor(0, i_NP, worker);

        int failed = worker.errIndex.load();
        if (failed >= 0) {                      // a worker task threw; ta_newP/ta_newC are incomplete
            Rcpp::stop("objective function failed for population member %d; "
                       "a compiled objective used with nthreads > 1 must not throw "
                       "or call the R API", failed + 1);
        }

        if (d_c > 0) {                        // JADE mean update (post-loop reduction, order-independent)
            double sumCR = 0.0, sumF = 0.0, sumF2 = 0.0;
            int nGood = 0;
            for (int i = 0; i < i_NP; i++) {
                if (accepted[i]) {
                    sumCR += rand_cross_used[i];
                    sumF  += rand_weight_used[i];
                    sumF2 += rand_weight_used[i] * rand_weight_used[i];
                    nGood++;
                }
            }
            if (nGood > 0) {
                double goodCR = sumCR / nGood;          // arithmetic mean of successful CR (standard JADE)
                meanCR = (1-d_c)*meanCR + d_c*goodCR;
                meanF  = (1-d_c)*meanF  + d_c*(sumF2/sumF);   // Lehmer mean of successful F, as in devol()
            }
        }

        if (i_bs_flag) {                        // best-of-parent-and-child selection, same semantics as devol()
            ta_popP.cols(0, i_NP-1)      = ta_oldP;     // reuse the caller's 2*NP scratch instead of
            ta_popC.rows(0, i_NP-1)      = ta_oldC;     // join_rows/join_cols allocating it per generation
            ta_popP.cols(i_NP, 2*i_NP-1) = ta_newP;
            ta_popC.rows(i_NP, 2*i_NP-1) = ta_newC;

            arma::uvec ord = arma::sort_index(ta_popC);
            for (int i = 0; i < i_NP; i++) {
                ta_newP.col(i) = ta_popP.col(ord[i]);
                ta_newC[i] = ta_popC[ord[i]];
            }
        }

        for (int i = 0; i < i_NP; i++) {        // best-tracking reduction (see function comment above)
            if (ta_newC[i] <= t_bestC) {
                t_bestC = ta_newC[i];
                t_bestP = ta_newP.col(i);
            }
        }

        ta_oldP.swap(ta_newP);                  // swap, not copy (see devol())
        ta_oldC.swap(ta_newC);

        if ( (i_trace > 0)  &&  ((i_iter % i_trace) == 0) ) {
            Rprintf("Iteration: %d/%d bestvalit: %f\n", i_iter, i_itermax, t_bestC);
        }

        if (std::abs(d_bestvalit[i_iter - 1] - t_bestC) < (d_reltol*(std::abs(d_bestvalit[i_iter - 1]) + d_reltol))) {
            i_iter_tol++;
        } else {
            i_iter_tol = 0;
        }
    } // end loop through generations

    d_pop = ta_oldP;
    i_iterations = i_iter;
    i_storepopcnt = popcnt;
    l_nfeval = static_cast<double>(ev.getNbEvals());
}
