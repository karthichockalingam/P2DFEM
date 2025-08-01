#include "mfem.hpp"

using namespace std;
using namespace mfem;

class ExchangeCurrentCoefficient: public Coefficient
{
   private:
      const GridFunction _surface_concentration;
      const GridFunction _electrolyte_concentration;

      GridFunctionCoefficient _sc;
      GridFunctionCoefficient _ec;
      TransformedCoefficient _jex;
   
   public:
      ExchangeCurrentCoefficient(
        const GridFunction & surface_concentration,
        const GridFunction & electrolyte_concentration):
        _surface_concentration(surface_concentration),
        _electrolyte_concentration(electrolyte_concentration),
        _sc(&_surface_concentration),
        _ec(&_electrolyte_concentration),
        _jex(&_sc, &_ec, [](real_t sc, real_t ec) { return sqrt( sc * ec * (1 - sc) ); }) {}

      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip)
      { return _jex.Eval(T, ip); }
};
