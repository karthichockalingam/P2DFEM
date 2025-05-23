#include "mfem.hpp"
using namespace mfem;

namespace constants {
    enum Method : int {
        SPM,
        SPMe,
        P2D
    };

    Method M = SPM;

    unsigned NPE = 0;
    unsigned NSEP = 0;
    unsigned NNE = 0;
    unsigned NX = 0;
    unsigned NPEPAR = 0;
    unsigned NNEPAR = 0;
    unsigned NPAR = 0;
    unsigned NR = 10;

    real_t LPE  = 0.;
    real_t LSEP = 0.;
    real_t LNE  = 0.;

    real_t AP = 1.;
    real_t AN = 1.;

    real_t DP = 1.;
    real_t DN = 1.;

    real_t CE0 = 1.;

    real_t I = 1.;

    void init_params(Method m) {
        M = m;

        switch (M)
        {
            case SPM:
                NPE = NSEP = NNE = 0;
                break;
            case SPMe:
                NPE = NNE = 2; NSEP = 1;
                break;
            case P2D:
                NPE = NSEP = NNE = 5;
                break;
        }
        NX = NPE + NSEP + NNE;

        switch (M)
        {
            case SPM:
                NPEPAR = NNEPAR = 1;
                LPE = LSEP = LNE = 1./3;
                break;
            case SPMe:
            case P2D:
                NPEPAR = NPE - 1;
                NNEPAR = NNE - 1;
                LPE  = 1. * NPE  / NX;
                LSEP = 1. * NSEP / NX;
                LNE  = 1. * NNE  / NX;
                break;
        }
        NPAR = NPEPAR + NNEPAR;
    }
}
