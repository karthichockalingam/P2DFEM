#include "mfem.hpp"
using namespace mfem;

class OpenCircuitPotentialCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _surface_concentration_gf;

      GridFunctionCoefficient _surface_concentration_gfc;

      const std::function<real_t(real_t)> _un;
      const std::function<real_t(real_t)> _up;

      const real_t & _scn;
      const real_t & _scp;

      TransformedCoefficient _ocp_ne_tc;
      TransformedCoefficient _ocp_pe_tc;

      PWConstCoefficient _ocp_pwcc;
      PWCoefficient _ocp_pwc;

   public:
      /// Default
      OpenCircuitPotentialCoefficient():
        _surface_concentration_gf(ParGridFunction()),
        _scn(0),
        _scp(0),
        _ocp_ne_tc(nullptr, _un),
        _ocp_pe_tc(nullptr, _up) {}

      /// SPM(e)
      OpenCircuitPotentialCoefficient(
        const std::function<real_t(real_t)> & un,
        const std::function<real_t(real_t)> & up,
        const real_t & scn,
        const real_t & scp):
        _surface_concentration_gf(ParGridFunction()),
        _un(un),
        _up(up),
        _scn(scn),
        _scp(scp),
        _ocp_ne_tc(nullptr, _un),
        _ocp_pe_tc(nullptr, _up),
        _ocp_pwcc(3) {}

      /// P2D
      OpenCircuitPotentialCoefficient(
        const std::function<real_t(real_t)> & un,
        const std::function<real_t(real_t)> & up,
        const ParGridFunction & sc):
        _surface_concentration_gf(sc),
        _surface_concentration_gfc(&_surface_concentration_gf),
        _un(un),
        _up(up),
        _scn(0),
        _scp(0),
        _ocp_ne_tc(&_surface_concentration_gfc, _un),
        _ocp_pe_tc(&_surface_concentration_gfc, _up),
        _ocp_pwc(Array<int>({NE, PE}), Array<Coefficient*>({static_cast<Coefficient*>(&_ocp_ne_tc), static_cast<Coefficient*>(&_ocp_pe_tc)})) {}

      /// SPM(e)
      virtual PWConstCoefficient & Eval()
      {
        _ocp_pwcc(NE) = _un(_scn);
        _ocp_pwcc(PE) = _up(_scp);
        return _ocp_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _ocp_pwc.Eval(T, ip);
      }

};
