#include "mfem.hpp"

using namespace std;
using namespace mfem;

class ExchangeCurrentCoefficient: public Coefficient
{
   private:
      const ParGridFunction _surface_concentration_gf;
      const ParGridFunction _electrolyte_concentration_gf;

      GridFunctionCoefficient _surface_concentration_gfc;
      GridFunctionCoefficient _electrolyte_concentration_gfc;

      ConstantCoefficient _jex_cc;
      TransformedCoefficient _jex_tc;
      Coefficient & _jex;

   public:
      /// SPM
      ExchangeCurrentCoefficient(
        const real_t & k,
        const real_t & sc,
        const real_t & ec):
        _jex_cc(k * sqrt( sc * ec * (1 - sc) )),
        _jex_tc(nullptr, nullptr, [](real_t, real_t) { return 0; }),
        _jex(_jex_cc) {}

      /// SPMe
      ExchangeCurrentCoefficient(
        const real_t & k,
        const real_t & sc,
        const ParGridFunction & ec):
        _electrolyte_concentration_gf(ec),
        _electrolyte_concentration_gfc(&_electrolyte_concentration_gf),
        _jex_tc(&_electrolyte_concentration_gfc, [=](real_t ec) { return k * sqrt( sc * ec * (1 - sc) ); }),
        _jex(_jex_tc) {}

      /// P2D
      ExchangeCurrentCoefficient(
        const real_t & k /* Probably needs to be a pwcoefficient */,
        const ParGridFunction & sc,
        const ParGridFunction & ec):
        _surface_concentration_gf(sc),
        _electrolyte_concentration_gf(ec),
        _surface_concentration_gfc(&_surface_concentration_gf),
        _electrolyte_concentration_gfc(&_electrolyte_concentration_gf),
        _jex_tc(&_surface_concentration_gfc, &_electrolyte_concentration_gfc, [=](real_t sc, real_t ec) { return k * sqrt( sc * ec * (1 - sc) ); }),
        _jex(_jex_tc) {}

      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _jex.Eval(T, ip);
      }

      virtual real_t Eval()
      {
        MFEM_ASSERT(&_jex == &_jex_cc, "ExchangeCurrentCoefficient does not wrap a ConstantCoefficient");
        return _jex_cc.constant;
      }
};
