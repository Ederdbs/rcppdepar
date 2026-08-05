
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2026  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#ifndef Rcpp_DE_evaluate_h_
#define Rcpp_DE_evaluate_h_

#include <Rcpp.h>
#include <atomic>
#include <cstdint>

namespace Rcpp {
    namespace DE {

        class EvalBase {
        public:
            EvalBase() : neval(0) {};
            virtual ~EvalBase() {}
            virtual double eval(NumericVector par) = 0;
            // uint64_t, not long: long is 32-bit on Win64 (LLP64) and 64-bit
            // elsewhere, which would overflow at ~2.1e9 evaluations on Windows
            // only and make the wrapped R type platform-dependent
            uint64_t getNbEvals() { return neval.load(std::memory_order_relaxed); }
        protected:
            // atomic: EvalCompiledRaw::evalRaw() is called concurrently from
            // RcppParallel worker threads when nthreads > 1
            std::atomic<uint64_t> neval;
        };

        class EvalStandard : public EvalBase {
        public:
            EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
            double eval(NumericVector par) {
                neval++;
                return defaultfun(par);
            }
        private:
            SEXP fcall, env;
            double defaultfun(NumericVector par) {        // essentialy same as the old evaluate
                Rcpp::Shield<SEXP> fn(::Rf_lang3(fcall, par, R_DotsSymbol));
                Rcpp::Shield<SEXP> sexp_fvec(::Rf_eval(fn, env));
                double f_result = REAL(sexp_fvec)[0];
                if (ISNAN(f_result))
                    Rcpp::stop("NaN value of objective function! \nPerhaps adjust the bounds.");
                return(f_result);
            }
        };

        typedef double (*funcPtr)(NumericVector);

        class EvalCompiled : public EvalBase {
        public:
            EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP env_) {
                funptr = *(xptr);
                env = env_;
            };
            EvalCompiled(SEXP xps, SEXP env_) {
                Rcpp::XPtr<funcPtr> xptr(xps);
                funptr = *(xptr);
                env = env_;
            };
            double eval(NumericVector par) {
                neval++;
                return funptr(par); //, env);
            }
        private:
            funcPtr funptr;
            SEXP env;
        };

        // Thread-safe objective: takes a raw double* + length instead of an
        // Rcpp::NumericVector, so it never touches R's (non-thread-safe)
        // SEXP allocator. This is the only objective type that may be called
        // from RcppParallel worker threads (see devol_parallel in devol.cpp).
        // Register such a function with putFunPtrInXPtrRaw()-style helpers.
        typedef double (*rawFuncPtr)(const double *, int);

        class EvalCompiledRaw : public EvalBase {
        public:
            EvalCompiledRaw(SEXP xps) {
                Rcpp::XPtr<rawFuncPtr> xptr(xps);
                funptr = *(xptr);
            }
            // main-thread entry point (keeps EvalBase interface usable when
            // nthreads == 1, e.g. from the existing serial devol() loop)
            double eval(NumericVector par) {
                return evalRaw(par.begin(), par.size());
            }
            // thread-safe entry point: no R API calls, safe to call
            // concurrently from RcppParallel worker threads
            double evalRaw(const double *par, int n) {
                neval.fetch_add(1, std::memory_order_relaxed);
                return funptr(par, n);
            }
        private:
            rawFuncPtr funptr;
        };

    }

}


#endif
