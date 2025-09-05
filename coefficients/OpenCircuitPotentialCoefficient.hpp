#include "mfem.hpp"
using namespace mfem;

class OpenCircuitPotentialCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _surface_concentration_gf;

      GridFunctionCoefficient _surface_concentration_gfc;

      TransformedCoefficient _ocp_ne_tc;
      ConstantCoefficient _ocp_sep_cc;
      TransformedCoefficient _ocp_pe_tc;

      Vector _ocp_vec;
      PWConstCoefficient _ocp_pwcc;
      PWCoefficient _ocp_pwc;
      Coefficient & _ocp;

   public:
      /// Default
      OpenCircuitPotentialCoefficient():
        _surface_concentration_gf(ParGridFunction()),
        _ocp_ne_tc(nullptr, [](real_t) { return 0; }),
        _ocp_pe_tc(nullptr, [](real_t) { return 0; }),
        _ocp(_ocp_pwcc) {}

      /// SPM(e)
      OpenCircuitPotentialCoefficient(
        const std::function<real_t(real_t)> & un,
        const std::function<real_t(real_t)> & up,
        const real_t & scn,
        const real_t & scp):
        _surface_concentration_gf(ParGridFunction()),
        _ocp_ne_tc(nullptr, [](real_t) { return 0; }),
        _ocp_pe_tc(nullptr, [](real_t) { return 0; }),
        _ocp_vec({un(scn), 0., up(scp)}),
        _ocp_pwcc(_ocp_vec),
        _ocp(_ocp_pwcc) {}

      /// P2D
      OpenCircuitPotentialCoefficient(
        const std::function<real_t(real_t)> & un,
        const std::function<real_t(real_t)> & up,
        const ParGridFunction & sc):
        _surface_concentration_gf(sc),
        _surface_concentration_gfc(&_surface_concentration_gf),
        _ocp_ne_tc(&_surface_concentration_gfc, un),
        _ocp_sep_cc(0),
        _ocp_pe_tc(&_surface_concentration_gfc, up),
        _ocp_pwc(Array<int>({NE, SEP, PE}), Array<Coefficient*>({static_cast<Coefficient*>(&_ocp_ne_tc), static_cast<Coefficient*>(&_ocp_sep_cc), static_cast<Coefficient*>(&_ocp_pe_tc)})),
        _ocp(_ocp_pwc) {}

      /// SPM(e)
      virtual PWConstCoefficient & Eval()
      {
        MFEM_ASSERT(&_ocp == &_ocp_pwcc, "OpenCircuitPotentialCoefficient does not wrap a PWConstCoefficient");
        return _ocp_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _ocp.Eval(T, ip);
      }

};
