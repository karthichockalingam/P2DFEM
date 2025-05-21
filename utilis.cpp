
#include "utilis.hpp"

const real_t m = 1.0;
const real_t cmax = 1.0;
const real_t T = 1.0;

double  FluxJExt(const double & ce, const double & cs)
{ 
   return m * pow(((cs - cmax) / (cs * ce)), 0.5);
}

double  FluxJ(const double & electrolyte_potential, const double & electrode_potential, const double & electrolyte_concentration, const double & electrode_surface_concentration, const Vector & x)
{ 
   double val = 2.0 * FluxJExt(electrolyte_concentration, electrode_surface_concentration) * (electrode_potential - electrolyte_potential) / (T);

   return std::sinh(val);
}

FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential, 
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration, 
    const double & electrode_surface_concentration, 
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
    const double & electrode_surface_concentration,   
    FuncT foo): 
    _electrolyte_potential(electrolyte_potential), 
    _electrode_potential(electrode_potential),
    _electrolyte_concentration(electrolyte_concentration),
    _electrode_surface_concentration(electrode_surface_concentration),
    TDGFunction( move(foo) ) {};

double 
FluxJGridFuncCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        return (GFunction)? GFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, transip) 
            : TDGFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), _electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, transip, GetTime() );
    }
/*
GridFuncFunctionCoefficient::GridFuncFunctionCoefficient (
    const GridFunction & electrolyte_concentration, 
    function<double(const double &)> foo): 
    _electrolyte_concentration(electrolyte_concentration),
    GFunction( move(foo) ) {};

double 
GridFuncFunctionCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        return GFunction(_electrolyte_concentration.GetValue(T, ip)); 
    }
*/

VectorGridFuncFunctionCoefficient::VectorGridFuncFunctionCoefficient(
    const GridFunction & electrolyte_concentration, 
    GradientGridFunctionCoefficient  & grad_electrolyte_concentration,
    ConstantCoefficient & kappa_D,
    function<real_t(const double &)> foo):
    VectorCoefficient(3),  
    _electrolyte_concentration(electrolyte_concentration),
    _grad_electrolyte_concentration(grad_electrolyte_concentration),
    _kappa_D(kappa_D),
    GFunction( move(foo) ) {};


void 
VectorGridFuncFunctionCoefficient::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        Vector  GradV;
        _grad_electrolyte_concentration.Eval(GradV, T, ip);
        // decision to be made: do we read in GradientGridFunctionCoefficient (like now) or create a
        // GradientGridFunctionCoefficient inside here each time?
        V = GradV;
        V *= _kappa_D.Eval(T, ip)/GFunction(_electrolyte_concentration.GetValue(T, ip)); 
    }
