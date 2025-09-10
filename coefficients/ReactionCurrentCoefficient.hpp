#include "mfem.hpp"
using namespace mfem;

class ReactionCurrentCoefficient: public Coefficient
{
   private:
      ExchangeCurrentCoefficient _jex;
      OverPotentialCoefficient _op;

      Vector _j_vec;
      PWConstCoefficient _j_pwcc;
      TransformedCoefficient _j_tc;
      Coefficient & _j;

   public:
      /// SPM(e)
      ReactionCurrentCoefficient():
        _j_vec({+ I / AN / LNE, 0., - I / AP / LPE}),
        _j_pwcc(_j_vec),
        _j_tc(nullptr, [](real_t) { return 0; }),
        _j(_j_pwcc) {}

      /// P2D
      ReactionCurrentCoefficient(
        const real_t & T,
        const ExchangeCurrentCoefficient & jex,
        const OverPotentialCoefficient & op):
        _jex(jex),
        _op(op),
        _j_tc(&_jex, &_op, [=](real_t jex, real_t op) { return 2 * jex * sinh (.5 * op / T); }),
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
