#include "mfem.hpp"

using namespace std;
using namespace mfem;

using Args = real_t(const real_t &, const real_t &, const real_t &, const real_t &);
using Func = function<Args>;


real_t  FluxJ(const real_t & electrolyte_potential, const real_t & electrode_potential, const real_t & electrolyte_concentration, const real_t & electrode_surface_concentration);
real_t  SPMeJFunc(const real_t & cs , const real_t & ce);

class FluxJGridFuncCoefficient : public Coefficient
 {
    const GridFunction & _electrolyte_potential;
    const GridFunction & _electrode_potential;
    const GridFunction & _electrolyte_concentration;
    const GridFunction & _electrode_surface_concentration;
    Func               GFunction;
 public:
    FluxJGridFuncCoefficient::FluxJGridFuncCoefficient(
    const GridFunction & electrolyte_potential,
    const GridFunction & electrode_potential,
    const GridFunction & electrolyte_concentration,
    const GridFunction & electrode_surface_concentration,
    Func foo):
    _electrolyte_potential(electrolyte_potential),
    _electrode_potential(electrode_potential),
    _electrolyte_concentration(electrolyte_concentration),
    _electrode_surface_concentration(electrode_surface_concentration),
    GFunction( move(foo) ) {};

    virtual   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };


 class ExchangeCurrentCoefficient: public Coefficient
 {
    const GridFunction & _surface_concentration;
    const GridFunction & _electrolyte_concentration;
 public:
    ExchangeCurrentCoefficient(
      const GridFunction & surface_concentration,
      const GridFunction & electrolyte_concentration):
      _surface_concentration(surface_concentration),
      _electrolyte_concentration(electrolyte_concentration) {};

    virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip);
 };
