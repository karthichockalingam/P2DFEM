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

    real_t Dp = 1e-14;
    real_t Dn = 3.9e-14;

    real_t r0 = 1e-6;
    real_t cpmax = 30555.;
    real_t cnmax = 51554.;
    real_t L = 3*80e-6; 
    real_t I_typ = 60.;

    real_t cell_area =  2.0527;
    real_t F = 96485.33289;

    real_t ts_p = F * cpmax * cell_area * L / I_typ;
    real_t ts_n = F * cnmax * cell_area * L / I_typ;

    real_t Dp_scale = pow(r0,2) / ts_p; // Diffusion coefficient scale for positive electrode (m^2/s).
    real_t Dn_scale = pow(r0,2) / ts_n; // Diffusion coefficient scale for negative electrode (m^2/s).

    real_t DP = Dp / Dp_scale; // Non-dimensional diffusion coefficient for positive electrode.
    real_t DN = Dn / Dn_scale; // Non-dimensional diffusion coefficient for negative electrode.

    real_t cp0 = cp0 / cpmax; 
    real_t cn0 = cn0 / cnmax;

    real_t CE0 = 1.;

    real_t I = 1.;

    void init_params(Method m, int order) {
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
                NPEPAR = NPE * order - 1;
                NNEPAR = NNE * order - 1;
                LPE  = 1. * NPE  / NX;
                LSEP = 1. * NSEP / NX;
                LNE  = 1. * NNE  / NX;
                break;
        }
        NPAR = NPEPAR + NNEPAR;
    }
}
