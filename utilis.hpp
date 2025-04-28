#include "mfem.hpp"

using namespace std;
using namespace mfem;

double  function3(const double & electrolyte_potential, const double & electrode_potential, const Vector & x);

class FluxJGridFuncCoefficient : public Coefficient
 {
    const GridFunction & _electrolyte_potential;
    const GridFunction & _electrode_potential;
    function<double(const double &, const double &, const Vector &)>                GFunction;
    function<double(const double &, const double &, const Vector &, const double)>  TDGFunction;
 public:
    FluxJGridFuncCoefficient(const GridFunction & electrolyte_potential, const GridFunction & electrode_potential
                       , function<double(const double &, const double &, const Vector &)> foo);
 
    FluxJGridFuncCoefficient(const GridFunction & electrolyte_potential, const GridFunction & electrode_potential
                       , function<double(const double &, const double &, const Vector &, const double)> foo);
 
    virtual   double Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };


 class FluxJExtGridFuncCoefficient : public Coefficient
 {
    const GridFunction & _electrolyte_concentration;
    const double & _electrode_surface_concentration;
    const double & _electrode_max_concentration;
    function<double(const double &, const double &, const double &, const Vector &)>                GFunction;
    function<double(const double &, const double &, const double &, const Vector &, const double)>  TDGFunction;
 public:
    FluxJExtGridFuncCoefficient(const GridFunction & electrolyte_concentration, const double & electrode_surface_concentration,
        const double & electrode_max_concentration, function<double(const double &, const double &, const double &, const Vector &)> foo);
 
    FluxJExtGridFuncCoefficient(GridFunction & electrolyte_concentration, const double & electrode_surface_concentration,
        const double & electrode_max_concentration, function<double(const double &, const double &, const double &, const Vector &, const double)> foo);
 
    virtual   double Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };
 