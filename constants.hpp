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
        PE,
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

    void init_params(Method m);
}
