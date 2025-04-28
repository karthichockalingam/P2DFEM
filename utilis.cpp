
#include "utilis.hpp"

double  function3(const double & electrolyte_potential, const double & electrode_potential, const Vector & x)
{ 
   double val = (electrode_potential - electrolyte_potential) / (2.0);

   return std::sinh(val);
}

double  function4(const double & ce, const double & cs, const double & cmax, const Vector & x)
{ 
   return pow(((cs - cmax) / (cs * ce)), 0.5);
}

FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(const GridFunction & electrolyte_potential, const GridFunction & electrode_potential
    , function<double(const double &, const double &, const Vector &)> foo)
    : _electrolyte_potential(electrolyte_potential), _electrode_potential(electrode_potential), 
    GFunction( move(foo) ) {};

FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(const GridFunction & electrolyte_potential, const GridFunction & electrode_potential
    , function<double(const double &, const double &, const Vector &, const double)> foo)
    : _electrolyte_potential(electrolyte_potential), _electrode_potential(electrode_potential),
    TDGFunction( move(foo) ) {};


double 
FluxJGridFuncCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        return (GFunction)? GFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), transip) 
            : TDGFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), transip, GetTime() );
    }

FluxJExtGridFuncCoefficient::FluxJExtGridFuncCoefficient(const GridFunction & electrolyte_concentration, const double & electrode_surface_concentration,
        const double & electrode_max_concentration, function<double(const double &, const double &, const double &, const Vector &)> foo)
        : _electrolyte_concentration(electrolyte_concentration), 
        _electrode_surface_concentration(electrode_surface_concentration),
        _electrode_max_concentration(electrode_max_concentration),
        GFunction( move(foo) ) {};
 
FluxJExtGridFuncCoefficient::FluxJExtGridFuncCoefficient(GridFunction & electrolyte_concentration, const double & electrode_surface_concentration,
        const double & electrode_max_concentration, function<double(const double &, const double &, const double &, const Vector &, const double)> foo)
        : _electrolyte_concentration(electrolyte_concentration), 
        _electrode_surface_concentration(electrode_surface_concentration),
        _electrode_max_concentration(electrode_max_concentration),
        TDGFunction( move(foo) ) {};
 
double 
FluxJExtGridFuncCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        return (GFunction)? GFunction(_electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, _electrode_max_concentration, transip) 
            : TDGFunction(_electrolyte_concentration.GetValue(T, ip), _electrode_surface_concentration, _electrode_max_concentration, transip, GetTime() );
    }