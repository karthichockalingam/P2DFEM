#include "mfem.hpp"
#include "parameters.hpp"

using namespace mfem;
using namespace LGM50;

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

    // Dimensional constants
    real_t F = 96485.33289; // Faraday constant, C/mol
    real_t R = 8.314; // Universal gas constant, J/(mol*K)
    real_t T_ref = 298.; // Reference temperature, K

    // Scalings
    real_t t0 = 1.0; // Time scale.
    real_t r0 = 1e-6;  // Length scale (particle)
    real_t L = positive_electrode_thickness + separator_thickness + negative_electrode_thickness; // Length scale (cell)

    real_t a0 = 1.0 / r0;

    real_t tp = F * cpmax * cell_area * L / I_typ; // Positive particle time scale.
    real_t tn = F * cnmax * cell_area * L / I_typ; // Negative particle time scale.

    real_t Dp_scale = r0 * r0 / t0; // Positive particle diffusion coefficient scale.  Units of m^2/s.
    real_t Dn_scale = r0 * r0 / t0; // Negative particle diffusion coefficient scale.  Units of m^2/s.
    // For some reason, JuBat calculates the following diffusion coefficient scales, but then later on multiplies
    // the mass matrix by tp / t0 (or tn / t0) which results in cancelling out tp (or tn) and replacing it with t0.
    //real_t Dp_scale = r0 * r0 / tp; // Positive particle diffusion coefficient scale.  Units of m^2/s.
    //real_t Dn_scale = r0 * r0 / tn; // Negative particle diffusion coefficient scale.  Units of m^2/s.

    real_t j_scale = I_typ / a0 / L / cell_area;

    // This works out to the same as j_scale above.
    //real_t j_scale_alt_p = r0 * cpmax * F / tp;
    //real_t j_scale_alt_n = r0 * cnmax * F / tn;

    real_t kp_scale = j_scale / cpmax / sqrt(ce0);
    real_t kn_scale = j_scale / cnmax / sqrt(ce0);

    real_t phi_scale = T_ref * R / F;

    // For scaling between particle time scale and cell time scale.  Required for scaling the flux j in the
    // SolidConcentration equation.
    real_t tp_scale = tp / t0;
    real_t tn_scale = tn / t0;

    // Scaled parameters
    real_t DP = Dp / Dp_scale; // Positive particle diffusion coefficient.
    real_t DN = Dn / Dn_scale; // Negative particle diffusion coefficient.

    real_t AP = Ap / a0; // Scaled positive electrode area.
    real_t AN = An / a0; // Scaled negative electrode area.

    real_t KP = kp_dim / kp_scale; // Scaled positive electrode reaction rate.
    real_t KN = kn_dim / kn_scale; // Scaled negative electrode reaction rate.

    real_t CP0 = cp0 / cpmax; // Scaled initial concentration, positive electrode.
    real_t CN0 = cn0 / cnmax; // Scaled initial concentration, negative electrode.

    real_t RP = rp / r0; // Scaled positive particle radius.
    real_t RN = rn / r0; // Scaled negative particle radius.

    // Extras to be properly defined later.
    real_t CE0 = 1.; // Scaled electrolyte concentration.
    real_t I = 1.; // Scaled current.
    real_t T = 1.0; // Scaled temperature.

    real_t UP(real_t ce) { return Up(ce); }
    real_t UN(real_t ce) { return Un(ce); }

    void init_params(Method m, int order) {
        M = m;

        switch (M)
        {
            case SPM:
                NPE = NSEP = NNE = 0;
                break;
            case SPMe:
            case P2D:
                NPE = NSEP = NNE = 5;
                break;
        }
        NX = NPE + NSEP + NNE;

        switch (M)
        {
            case SPM:
            case SPMe:
                NPEPAR = NNEPAR = 1;
                LPE = positive_electrode_thickness / L;
                LSEP = separator_thickness / L;
                LNE = negative_electrode_thickness / L;
                break;
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
