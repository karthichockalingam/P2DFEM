#include "mfem.hpp"

using namespace mfem;

namespace constants {
    enum PotentialBlock : int {
        EPP,
        SPP,
    };

    enum ConcentrationBlock : int {
        ECC,
        SCC
    };

    enum Block : int {
        P = 0,
        C = 2
    };

    enum XBlock : int {
        EP = P + EPP,
        SP = P + SPP,
        EC = C + ECC,
        SC = C + SCC
    };

    enum Region : int {
        UNKNOWN,
        NE,
        SEP,
        PE
    };

    enum Model : int {
        SPM,
        SPMe,
        P2D
    };

    extern const Model M;

    extern const unsigned NNE;
    extern const unsigned NSEP;
    extern const unsigned NPE;
    extern const unsigned NX;
    extern const unsigned NR;

    extern const unsigned NNEPAR;
    extern const unsigned NPEPAR;
    extern const unsigned NPAR;

    extern const unsigned NMACROP;
    extern const unsigned NMACROC;
    extern const unsigned NMACRO;
    extern const unsigned NEQS;

    extern const real_t LNE;
    extern const real_t LSEP;
    extern const real_t LPE;

    extern const real_t AN;
    extern const real_t AP;

    extern const real_t DN;
    extern const real_t DP;

    extern const real_t KN;
    extern const real_t KP;

    extern const real_t RN;
    extern const real_t RP;

    extern const real_t CN0;
    extern const real_t CP0;

    extern const real_t CE0;

    extern const real_t I;
    extern const real_t T;

    extern const real_t phi_scale;
    extern const real_t tn_scale;
    extern const real_t tp_scale;

    real_t UN(real_t);
    real_t UP(real_t);

    void init_params(std::string m, int order);
}
