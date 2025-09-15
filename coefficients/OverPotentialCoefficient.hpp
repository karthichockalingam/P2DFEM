#include "mfem.hpp"
using namespace mfem;

class OverPotentialCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _solid_potential_gf;
      const ParGridFunction & _electrolyte_potential_gf;

      GridFunctionCoefficient _solid_potential_gfc;
      GridFunctionCoefficient _electrolyte_potential_gfc;

      ExchangeCurrentCoefficient * _jex;
      OpenCircuitPotentialCoefficient * _ocp;

      SumCoefficient _dp_sc;

      PWConstCoefficient _op_pwcc;
      SumCoefficient _op_sc;

   public:
      /// Default
      OverPotentialCoefficient():
        _solid_potential_gf(ParGridFunction()),
        _electrolyte_potential_gf(ParGridFunction()),
        _dp_sc(0, _solid_potential_gfc),
        _op_sc(0, _solid_potential_gfc) {}

      /// SPM(e)
      OverPotentialCoefficient(
        const real_t & T,
        ExchangeCurrentCoefficient & jex):
        _solid_potential_gf(ParGridFunction()),
        _electrolyte_potential_gf(ParGridFunction()),
        _dp_sc(0, _solid_potential_gfc),
        _op_sc(0, _solid_potential_gfc),
        _jex(&jex),
        _op_pwcc(3) {}

      /// P2D
      OverPotentialCoefficient(
        const ParGridFunction & sp,
        const ParGridFunction & ep,
        OpenCircuitPotentialCoefficient & ocp):
        _solid_potential_gf(sp),
        _electrolyte_potential_gf(ep),
        _solid_potential_gfc(&_solid_potential_gf),
        _electrolyte_potential_gfc(&_electrolyte_potential_gf),
        _ocp(&ocp),
        _dp_sc(_solid_potential_gfc, _electrolyte_potential_gfc, 1, -1),
        _op_sc(_dp_sc, *_ocp, 1, -1) {}

      /// SPM(e)
      virtual PWConstCoefficient & Eval()
      {
        _op_pwcc(NE) = 2 * T * asinh(+ I / AN / LNE / 2.0 / _jex->Eval()(NE));
        _op_pwcc(PE) = 2 * T * asinh(- I / AP / LPE / 2.0 / _jex->Eval()(PE));
        return _op_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _op_sc.Eval(T, ip);
      }

};
