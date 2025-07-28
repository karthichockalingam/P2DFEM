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

    // Dimensional constants
    real_t F = 96485.33289; // Faraday constant, C/mol
    
    

    // Electrochemical parameters
    real_t cpmax = 63104; // Maximum concentration, positive electrode, mol/m^3
    real_t cnmax = 33133; // Maximum concentration, negative electrode, mol/m

    real_t cell_length = 1.58;
    real_t cell_width = 6.5e-2;
    real_t cell_no_layers = 1;
    real_t cell_area = cell_width * cell_length * cell_no_layers;

    real_t I_typ = 5.0;    // Or I1C in Jubat.

    real_t eps_p = 0.335; // Porosity???
    real_t eps_p_fi = 0.0;
    real_t eps_p_s = 1 - eps_p - eps_p_fi;

    real_t eps_n = 0.25; // Porosity???
    real_t eps_n_fi = 0.0;
    real_t eps_n_s = 1 - eps_n - eps_n_fi;

    real_t rs_p = 5.22e-6; // Positive particle radius (m).
    real_t rs_n = 5.86e-6; // Negative particle radius (m).

    real_t positive_electrode_thickness = 75.6e-6;
    real_t separator_thickness = 12e-6;
    real_t negative_electrode_thickness = 85.2e-6;

    real_t kp_dim = 3.42e-6;
    real_t kn_dim = 6.48e-7;

    real_t ce0 = 1000.0;


    // Scalings
    real_t r0 = 1e-6;  // Length scale (particle)
    real_t L = 1e-6;  // Length scale (cell)
    
    real_t a0 = 1.0 / r0; 

    real_t tp = F * cpmax * cell_area * L / I_typ;
    real_t tn = F * cnmax * cell_area * L / I_typ;

    real_t Dp_scale = r0 * r0 / tp; // Positive particle diffusion coefficient scale.
    real_t Dn_scale = r0 * r0 / tn; // Negative particle diffusion coefficient scale.

    real_t j_scale = I_typ / a0 / L / cell_area;

    real_t kp_scale = j_scale / cpmax / sqrt(ce0);
    real_t kn_scale = j_scale / cnmax / sqrt(ce0);


    // Dimensional parameters
    real_t Dp = 4.0e-15; // Diffusion coefficient, positive electrode
    real_t Dn = 3.3e-14; // Diffusion coefficient, negative electrode

    real_t Ap = 3 * eps_p_s / rs_p; // Positive electrode area (m^2).
    real_t An = 3 * eps_n_s / rs_n; // Negative electrode area (m^2).


    // Scaled parameters
    real_t DP = Dp / Dp_scale; // Positive particle diffusion coefficient.
    real_t DN = Dn / Dn_scale; // Negative particle diffusion coefficient.

    real_t AP = Ap / a0; // Scaled positive electrode area.
    real_t AN = An / a0; // Scaled negative electrode area.

    real_t KP = kp_dim / kp_scale; // Scaled positive electrode reaction rate.
    real_t KN = kn_dim / kn_scale; // Scaled negative electrode reaction rate.

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
                //LPE = LSEP = LNE = 1./3;
                LPE = positive_electrode_thickness / L;
                LSEP = separator_thickness / L;
                LNE = negative_electrode_thickness / L;
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
