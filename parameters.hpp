#pragma once
#include "mfem.hpp"

using namespace mfem;

namespace LGM50 {
    // Positive electrode.
    const real_t cp0 = 17038.; // Initial concentration, positive electrode, mol/m^3
    const real_t cpmax = 63104; // Maximum concentration, positive electrode, mol/m^3

    const real_t rp = 5.22e-6; // Positive particle radius (m).
    const real_t positive_electrode_thickness = 75.6e-6;

    const real_t eps_p = 0.335; // Porosity???
    const real_t eps_p_fi = 0.0;
    const real_t eps_p_s = 1 - eps_p - eps_p_fi;

    const real_t kp_dim = 3.42e-6;
    const real_t Dp = 4.0e-15; // Diffusion coefficient, positive electrode
    const real_t Ap = 3 * eps_p_s / rp; // Positive electrode area (m^2).

    real_t Up(real_t cs) // Open circuit potential (V)
    {
        return -0.8090*cs + 4.4875 - 0.0428*tanh(18.5138*(cs - 0.5542))
               -17.7326*tanh(15.7890*(cs - 0.3117)) + 17.5842*tanh(15.9308*(cs - 0.3120));
    }

    // Negative electrode.
    const real_t cn0 = 29866.; // Initial concentration, negative electrode, mol/m^3
    const real_t cnmax = 33133; // Maximum concentration, negative electrode, mol/m

    const real_t rn = 5.86e-6; // Negative particle radius (m).
    const real_t negative_electrode_thickness = 85.2e-6;

    const real_t eps_n = 0.25; // Porosity???
    const real_t eps_n_fi = 0.0;
    const real_t eps_n_s = 1 - eps_n - eps_n_fi;

    const real_t kn_dim = 6.48e-7;
    const real_t Dn = 3.3e-14; // Diffusion coefficient, negative electrode
    const real_t An = 3 * eps_n_s / rn; // Negative electrode area (m^2).

    real_t Un(real_t cs) // Open circuit potential (V)
    {
        return 1.97938*exp(-39.3631*cs) + 0.2482 - 0.0909*tanh(29.8538*(cs - 0.1234))
               -0.04478*tanh(14.9159*(cs - 0.2769)) - 0.0205*tanh(30.4444*(cs - 0.6103));
    }

    // Separator.
    const real_t separator_thickness = 12e-6;

    // Electrolyte.
    const real_t ce0 = 1000.0;

    // Cell parameters.
    const real_t cell_length = 1.58;
    const real_t cell_width = 6.5e-2;
    const real_t cell_no_layers = 1;
    const real_t cell_area = cell_width * cell_length * cell_no_layers;

    const real_t I_typ = 5.0;    // Or I1C in Jubat.
}
