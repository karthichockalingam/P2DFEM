#include "mfem.hpp"

using namespace std;
using namespace mfem;

using Args = real_t(const real_t &, const real_t &, const real_t &, const real_t &, const Vector &);
using Func = function<Args>;
using ArgsT = real_t(const real_t &, const real_t &, const real_t &, const real_t &, const Vector &, const real_t);
using FuncT = function<ArgsT>;

real_t  FluxJExt(const real_t & ce, const real_t & cs);
real_t  FluxJ(const real_t & electrolyte_potential, const real_t & electrode_potential, const real_t & electrolyte_concentration, const real_t & electrode_surface_concentration, const Vector & x);
real_t  SPMeJFunc(const real_t & cs , const real_t & ce , const Vector & x);

class FluxJGridFuncCoefficient : public Coefficient
 {
    const GridFunction & _electrolyte_potential;
    const GridFunction & _electrode_potential;
    const GridFunction & _electrolyte_concentration;
    const real_t & _electrode_surface_concentration;
    Func               GFunction;
    FuncT              TDGFunction;
 public:
    FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential,
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const real_t & _electrode_surface_concentration,
    Func foo);

    FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential,
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const real_t & _electrode_surface_concentration,
    FuncT foo);

    virtual   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };


 class ExternalCurrentCoefficient: public Coefficient
 {
    const GridFunction & _surface_concentraion;
    const GridFunction & _electrolyte_concentration;
 public:
    ExternalCurrentCoefficient(
      const GridFunction & surface_concentraion,
      const GridFunction & electrolyte_concentration):
      _surface_concentraion(surface_concentraion),
      _electrolyte_concentration(electrolyte_concentration) {};

    virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };
