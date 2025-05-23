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

    extern Method M;

    extern unsigned NPE;
    extern unsigned NSEP;
    extern unsigned NNE;
    extern unsigned NX;
    extern unsigned NPAR;
    extern unsigned NR;

    extern real_t LPE;
    extern real_t LSEP;
    extern real_t LNE;

    void init_params();
}
