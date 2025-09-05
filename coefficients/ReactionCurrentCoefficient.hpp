#include "mfem.hpp"
using namespace mfem;

class ReactionCurrentCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _solid_potential_gf;
      const ParGridFunction & _electrolyte_potential_gf;

      GridFunctionCoefficient _solid_potential_gfc;
      GridFunctionCoefficient _electrolyte_potential_gfc;

      ExchangeCurrentCoefficient _jex;
      OpenCircuitPotentialCoefficient _ocp;

      SumCoefficient _dp_sc;
      SumCoefficient _op_sc;

      Vector _j_vec;
      PWConstCoefficient _j_pwcc;
      TransformedCoefficient _j_tc;
      Coefficient & _j;

   public:
      /// SPM(e)
      ReactionCurrentCoefficient():
        _solid_potential_gf(ParGridFunction()),
        _electrolyte_potential_gf(ParGridFunction()),
        _dp_sc(0, _solid_potential_gfc),
        _op_sc(0, _solid_potential_gfc),
        _j_vec({+ I / AN / LNE, 0., - I / AP / LPE}),
        _j_pwcc(_j_vec),
        _j_tc(nullptr, [](real_t) { return 0; }),
        _j(_j_pwcc) {}

      /// P2D
      ReactionCurrentCoefficient(
        const real_t & T,
        const ParGridFunction & sp,
        const ParGridFunction & ep,
        const ExchangeCurrentCoefficient & jex,
        const OpenCircuitPotentialCoefficient & ocp):
        _solid_potential_gf(sp),
        _electrolyte_potential_gf(ep),
        _solid_potential_gfc(&_solid_potential_gf),
        _electrolyte_potential_gfc(&_electrolyte_potential_gf),
        _jex(jex),
        _ocp(ocp),
        _dp_sc(_solid_potential_gfc, _solid_potential_gfc, 1, -1),
        _op_sc(_dp_sc, _ocp, 1, -1),
        _j_tc(&_jex, &_op_sc, [=](real_t jex, real_t op) { return 2 * jex * sinh (.5 * op / T); }),
        _j(_j_tc) {}

      /// SPM(e)
      virtual PWConstCoefficient Eval()
      {
        MFEM_ASSERT(&_j == &_j_pwcc, "ReactionCurrentCoefficient does not wrap a PWConstCoefficient");
        return _j_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _j.Eval(T, ip);
      }

};
