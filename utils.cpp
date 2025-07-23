
#include "utils.hpp"

const real_t m = 1.0;
const real_t cmax = 1.0;
const real_t T = 1.0;

real_t SPMeJFunc(const real_t & cs, const real_t & ce)
{
   return m * sqrt(ce) * sqrt(cs) * sqrt(cs - cmax);
}

real_t  FluxJ(const real_t & electrolyte_potential, const real_t & electrode_potential, const real_t & electrolyte_concentration, const real_t & electrode_surface_concentration)
{
   real_t val = 2.0 * SPMeJFunc(electrolyte_concentration, electrode_surface_concentration) * (electrode_potential - electrolyte_potential) / (T);

   return std::sinh(val);
}

real_t
FluxJGridFuncCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return GFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration.GetValue(T, ip));
}

real_t
ExchangeCurrentCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return SPMeJFunc(_surface_concentration.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip));
}
