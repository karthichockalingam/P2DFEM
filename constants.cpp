#include "mfem.hpp"
#include "parameters.hpp"

using namespace mfem;
using namespace LGM50;

namespace constants {
    enum Model : int {
        SPM,
        SPMe,
        P2D
    };

    Model M = SPM;

    unsigned NNE = 0;
    unsigned NSEP = 0;
    unsigned NPE = 0;
    unsigned NX = 0;
    unsigned NR = 10;

    unsigned NNEPAR = 0;
    unsigned NPEPAR = 0;
    unsigned NPAR = 0;

    unsigned NMACROP = 2;
    unsigned NMACROC = 1;
    unsigned NMACRO = NMACROP + NMACROC;
    unsigned NEQS = 0;

    // Dimensional constants
    real_t F = 96485.33289; // Faraday constant, C/mol
    real_t R = 8.314; // Universal gas constant, J/(mol*K)
    real_t T_ref = 298.; // Reference temperature, K

    // Scalings
    real_t t0 = 1.0; // Time scale.
    real_t r0 = 1e-6;  // Length scale (particle)
    real_t L = negative_electrode_thickness + separator_thickness + positive_electrode_thickness; // Length scale (cell)

    real_t LNE  = negative_electrode_thickness / L;
    real_t LSEP = separator_thickness          / L;
    real_t LPE  = positive_electrode_thickness / L;

    real_t a0 = 1.0 / r0;

    real_t tn = F * cnmax * cell_area * L / I_typ; // Negative particle time scale.
    real_t tp = F * cpmax * cell_area * L / I_typ; // Positive particle time scale.

    real_t Dn_scale = r0 * r0 / t0; // Negative particle diffusion coefficient scale.  Units of m^2/s.
    real_t Dp_scale = r0 * r0 / t0; // Positive particle diffusion coefficient scale.  Units of m^2/s.
    // For some reason, JuBat calculates the following diffusion coefficient scales, but then later on multiplies
    // the mass matrix by tp / t0 (or tn / t0) which results in cancelling out tp (or tn) and replacing it with t0.
    //real_t Dn_scale = r0 * r0 / tn; // Negative particle diffusion coefficient scale.  Units of m^2/s.
    //real_t Dp_scale = r0 * r0 / tp; // Positive particle diffusion coefficient scale.  Units of m^2/s.

    real_t j_scale = I_typ / a0 / L / cell_area;

    // This works out to the same as j_scale above.
    //real_t j_scale_alt_n = r0 * cnmax * F / tn;
    //real_t j_scale_alt_p = r0 * cpmax * F / tp;

    real_t kn_scale = j_scale / cnmax / sqrt(ce0);
    real_t kp_scale = j_scale / cpmax / sqrt(ce0);

    real_t phi_scale = T_ref * R / F;

    // For scaling between particle time scale and cell time scale.  Required for scaling the flux j in the
    // SolidConcentration equation.
    real_t tn_scale = tn / t0;
    real_t tp_scale = tp / t0;

    // Scaled parameters
    real_t DN = Dn / Dn_scale; // Diffusion coefficient of each Negative particle
    real_t DP = Dp / Dp_scale; // Diffusion coefficient of each Negative particle

    real_t AN = An / a0; // scaled surface Area of each Negative particle // later SAN!
    real_t AP = Ap / a0; // scaled surface Area of each Positive particle

    real_t KN = kn_dim / kn_scale; // scaled reaction rate of each Negative particle
    real_t KP = kp_dim / kp_scale; // scaled reaction rate of each Positive particle

    real_t CN0 = cn0 / cnmax;  // scaled initial Concentration of Negative particle
    real_t CP0 = cp0 / cpmax; // scaled initial Concentration of Positive particle

    real_t RN = rn / r0; // scaled Radius of Negative particle
    real_t RP = rp / r0; // scaled Radius of Positive particle

    // Extras to be properly defined later.
    real_t CE0 = 1.; // scaled initial Concentration of Electrolyte
    real_t I = 1.; // scaled external current
    real_t T = 1.0; // scaled Temperature.

    real_t UN(real_t ce) { return Un(ce); }
    real_t UP(real_t ce) { return Up(ce); }

    void init_params(std::string m, int order) {
        std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c){ return std::tolower(c); });

        if (m == "spm")
            M = SPM;
        else if (m == "spme")
            M = SPMe;
        else if (m == "p2d" || m == "dfn")
            M = P2D;
        else
            MFEM_ASSERT(false, "Unrecognised model.");

        switch (M)
        {
            case SPM:
                NNE = NSEP = NPE = 0;
                break;
            case SPMe:
            case P2D:
                NNE = NSEP = NPE = 5;
                break;
        }
        NX = NNE + NSEP + NPE;

        switch (M)
        {
            case SPM:
            case SPMe:
                NNEPAR = NPEPAR = 1;
                break;
            case P2D:
                NNEPAR = NNE * order - 1;
                NPEPAR = NPE * order - 1;
                break;
        }
        NPAR = NNEPAR + NPEPAR;
        NEQS = NMACRO + NPAR;
    }
}
