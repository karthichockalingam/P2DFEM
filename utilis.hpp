#include "mfem.hpp"

using namespace std;
using namespace mfem;

using Args = double(const double &, const double &, const double &, const double &, const Vector &);
using Func = function<Args>;
using ArgsT = double(const double &, const double &, const double &, const double &, const Vector &, const double);
using FuncT = function<ArgsT>;

double  FluxJExt(const double & ce, const double & cs);
double  FluxJ(const double & electrolyte_potential, const double & electrode_potential, const double & electrolyte_concentration, const double & electrode_surface_concentration, const Vector & x);

class FluxJGridFuncCoefficient : public Coefficient
 {
    const GridFunction & _electrolyte_potential;
    const GridFunction & _electrode_potential;
    const GridFunction & _electrolyte_concentration;
    const double & _electrode_surface_concentration;
    Func               GFunction;
    FuncT              TDGFunction;
 public:
    FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential, 
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const double & _electrode_surface_concentration, 
    Func foo);
 
    FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential, 
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration, 
    const double & _electrode_surface_concentration, 
    FuncT foo);
 
    virtual   double Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };



 