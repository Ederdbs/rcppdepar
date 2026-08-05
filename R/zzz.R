.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste("Now loading:",
                                "  RcppDEpar: parallel C++ Implementation of Differential Evolution Optimisation",
                                "  Author: Eder David Borges da Silva",
                                "  Based on:",
                                "  RcppDE: C++ Implementation of Differential Evolution Optimisation",
                                "  Author: Dirk Eddelbuettel",
                                "  DEoptim (version 2.0-7): Differential Evolution algorithm in R",
                                "  Authors: David Ardia, Katharine Mullen, Brian Peterson and Joshua Ulrich",
                                sep="\n"))
}
