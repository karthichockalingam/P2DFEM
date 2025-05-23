#include "constants.hpp"

namespace constants {
    Method M = SPM;

    unsigned NPE;
    unsigned NSEP;
    unsigned NNE;
    unsigned NX;
    unsigned NPAR;
    unsigned NR = 10;

    real_t LPE;
    real_t LSEP;
    real_t LNE;

    void init_params() {
        switch (M)
        {
            case SPM:
            case SPMe:
                NPE = NSEP = NNE = NX = 0;
                NPAR = 2;
                break;
            case P2D:
                NPE = NSEP = NNE = 5;
                NX = NPE + NSEP + NNE;
                NPAR = (NPE - 1) + (NNE - 1);
                break;
        }
        LPE  = 1. * NPE  / NX;
        LSEP = 1. * NSEP / NX;
        LNE  = 1. * NNE  / NX;
    }
}
