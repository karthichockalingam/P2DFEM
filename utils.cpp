
#include "utils.hpp"

const real_t m = 1.0;
const real_t cmax = 1.0;
const real_t T = 1.0;

real_t  FluxJExt(const real_t & ce, const real_t & cs)
{
   return m * pow(((cs - cmax) / (cs * ce)), 0.5);
}

real_t  FluxJ(const real_t & electrolyte_potential, const real_t & electrode_potential, const real_t & electrolyte_concentration, const real_t & electrode_surface_concentration, const Vector & x)
{
   real_t val = 2.0 * FluxJExt(electrolyte_concentration, electrode_surface_concentration) * (electrode_potential - electrolyte_potential) / (T);

   return std::sinh(val);
}

real_t SPMeJFunc(const real_t & cs, const real_t & ce)
{
   return m * sqrt(ce) * sqrt(cs) * sqrt(cs - cmax);
}

FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential,
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const real_t & electrode_surface_concentration,
    Func foo):
    _electrolyte_potential(electrolyte_potential),
    _electrode_potential(electrode_potential),
    _electrolyte_concentration(electrolyte_concentration),
    _electrode_surface_concentration(electrode_surface_concentration),
    GFunction( move(foo) ) {};

FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential,
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const real_t & electrode_surface_concentration,
    FuncT foo):
    _electrolyte_potential(electrolyte_potential),
    _electrode_potential(electrode_potential),
    _electrolyte_concentration(electrolyte_concentration),
    _electrode_surface_concentration(electrode_surface_concentration),
    TDGFunction( move(foo) ) {};

real_t
FluxJGridFuncCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        real_t x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        return (GFunction)? GFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, transip)
            : TDGFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, transip, GetTime() );
    }

real_t
ExternalCurrentCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        Vector transip(3);
        T.Transform(ip, transip);
        return SPMeJFunc(_surface_concentraion.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip));
    }

