#include "mfem.hpp"
#include "cells/LGM50.hpp"

using namespace mfem;
using namespace LGM50;

namespace constants {
    bool SPM = false;
    bool SPMe = false;
    bool P2D = false;

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
    real_t r0 = 1e-6;  // Length scale (particle), m
    real_t L = negative_electrode_thickness + separator_thickness + positive_electrode_thickness; // Length scale (cell), m

    real_t LNE  = negative_electrode_thickness / L;
    real_t LSEP = separator_thickness          / L;
    real_t LPE  = positive_electrode_thickness / L;

    real_t a0 = 1.0 / r0; // specific surface area scale, m^{-1}

    real_t tn = F * cnmax * cell_area * L / I_typ; // Negative particle time scale.
    real_t tp = F * cpmax * cell_area * L / I_typ; // Positive particle time scale.

    real_t te = F * ce0 * cell_area * L / I_typ; // Electrolyte "particle" time scale.

    real_t te_scale = te / t0;

    real_t Dn_scale = r0 * r0 / t0; // Negative particle diffusion coefficient scale.  Units of m^2/s.
    real_t Dp_scale = r0 * r0 / t0; // Positive particle diffusion coefficient scale.  Units of m^2/s.
    // For some reason, JuBat calculates the following diffusion coefficient scales, but then later on multiplies
    // the mass matrix by tp / t0 (or tn / t0) which results in cancelling out tp (or tn) and replacing it with t0.
    //real_t Dn_scale = r0 * r0 / tn; // Negative particle diffusion coefficient scale.  Units of m^2/s.
    //real_t Dp_scale = r0 * r0 / tp; // Positive particle diffusion coefficient scale.  Units of m^2/s.

    real_t De_scale = L * L / te;

    // Transport efficiency (inverse MacMullin number). This is B(x) in Planella, and is absorbed into the
    // definition of kappa_ne/kappa_pe/kappa_sp in JuBat.
    real_t BPE = pow(eps_p, brugg);
    real_t BNE = pow(eps_n, brugg);
    real_t BSEP = pow(eps_s, brugg);

    real_t j_scale = I_typ / a0 / L / cell_area; //volumetric current density scale, A/m^2

    // This works out to the same as j_scale above.
    //real_t j_scale_alt_n = r0 * cnmax * F / tn;
    //real_t j_scale_alt_p = r0 * cpmax * F / tp;

    real_t kn_scale = j_scale / cnmax / sqrt(ce0);
    real_t kp_scale = j_scale / cpmax / sqrt(ce0);

    real_t phi_scale = T_ref * R / F; //thermal voltage scale, V

    // For scaling between particle time scale and cell time scale.  Required for scaling the flux j in the
    // SolidConcentration equation.
    real_t tn_scale = tn / t0;
    real_t tp_scale = tp / t0;

    real_t ce_scale = ce0;

    real_t sig_scale = L * I_typ / (cell_area * phi_scale); // Electrode conductivity scale.
    real_t kappa_scale = L * I_typ / (cell_area * phi_scale); // Electrolyte conductivity scale.

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

    real_t SIGP = sig_p / sig_scale; // Scaled positive electrode conductivity.
    real_t SIGN = sig_n / sig_scale; // Scaled negative electrode conductivity.

    // Extras to be properly defined later.
    real_t CE0 = ce0 / ce_scale; // Scaled electrolyte concentration.
    real_t I = 1.; // scaled external current
    real_t T = 1.0; // scaled Temperature.

    real_t EPS_P = eps_p;
    real_t EPS_N = eps_n;
    real_t EPS_S = eps_s;

    real_t TPLUS = tplus;

    real_t UN(real_t ce) { return Un(ce) / phi_scale; }
    real_t UP(real_t ce) { return Up(ce) / phi_scale; }
    real_t DE(real_t ce) { return De(ce * ce_scale) / De_scale; }
    real_t Kappa(real_t x) { return kappa(x * ce_scale) / kappa_scale; }

    real_t KS = Kappa(CE0);// / kappa_scale; // Scaled electrolyte conductivity.

    void init_params(std::string m, int order)
    {
        std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c){ return std::tolower(c); });

        if (m == "spm")
            SPM = true;
        else if (m == "spme")
            SPMe = true;
        else if (m == "p2d" || m == "dfn")
            P2D = true;
        else
            MFEM_ASSERT(false, "Unrecognised model.");

        if (SPM)
            NNE = NSEP = NPE = 0;
        else if (SPMe || P2D)
            NNE = NSEP = NPE = 10;

        NX = NNE + NSEP + NPE;

        if (SPM || SPMe)
            NNEPAR = NPEPAR = 1;
        else if (P2D)
        {
            NNEPAR = NNE * order + 1;
            NPEPAR = NPE * order + 1;
        }
        NPAR = NNEPAR + NPEPAR;
        NEQS = NMACRO + NPAR;
    }
}
