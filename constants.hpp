#include "mfem.hpp"

using namespace mfem;

namespace constants {
    enum Block : int {
        EP,
        EC,
        SP,
        SC
    };

    enum Region : int {
        PE = 1,
        SEP,
        NE
    };

    enum Method : int {
        SPM,
        SPMe,
        P2D
    };

    extern const Method M;

    extern const unsigned NPE;
    extern const unsigned NSEP;
    extern const unsigned NNE;
    extern const unsigned NX;
    extern const unsigned NPEPAR;
    extern const unsigned NNEPAR;
    extern const unsigned NPAR;
    extern const unsigned NR;

    extern const real_t LPE;
    extern const real_t LSEP;
    extern const real_t LNE;

    extern const real_t AP;
    extern const real_t AN;

    extern const real_t DP;
    extern const real_t DN;

    extern const real_t KP;
    extern const real_t KN;

    extern const real_t RP;
    extern const real_t RN;

    extern const real_t CP0;
    extern const real_t CN0;

    extern const real_t CE0;

    extern const real_t I;
    extern const real_t T;

    extern const real_t phi_scale;
    extern const real_t tp_scale;
    extern const real_t tn_scale;

    real_t UP(real_t);
    real_t UN(real_t);

    void init_params(Method m, int order);
}
