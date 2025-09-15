#include "mfem.hpp"
using namespace mfem;

class ExchangeCurrentCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _surface_concentration_gf;
      const ParGridFunction & _electrolyte_concentration_gf;

      GridFunctionCoefficient _surface_concentration_gfc;
      GridFunctionCoefficient _electrolyte_concentration_gfc;

      const real_t & _scn;
      const real_t & _scp;

      TransformedCoefficient _jex_ne_tc;
      TransformedCoefficient _jex_pe_tc;

      Vector _jex_vec;
      PWConstCoefficient _jex_pwcc;
      PWCoefficient _jex_pwc;

   public:
      /// Default
      ExchangeCurrentCoefficient():
        _surface_concentration_gf(ParGridFunction()),
        _electrolyte_concentration_gf(ParGridFunction()),
        _scn(0),
        _scp(0),
        _jex_ne_tc(nullptr, [](real_t) { return 0; }),
        _jex_pe_tc(nullptr, [](real_t) { return 0; }) {}

      /// SPM
      ExchangeCurrentCoefficient(
        const real_t & kn,
        const real_t & kp,
        const real_t & scn,
        const real_t & scp,
        const real_t & ec):
        _surface_concentration_gf(ParGridFunction()),
        _electrolyte_concentration_gf(ParGridFunction()),
        _scn(scn),
        _scp(scp),
        _jex_ne_tc(nullptr, [](real_t) { return 0; }),
        _jex_pe_tc(nullptr, [](real_t) { return 0; }),
        _jex_vec({kn * sqrt( ec ), 0., kp * sqrt( ec )}),
        _jex_pwcc(3) {}

      /// SPMe
      ExchangeCurrentCoefficient(
        const real_t & kn,
        const real_t & kp,
        const real_t & scn,
        const real_t & scp,
        const ParGridFunction & ec):
        _surface_concentration_gf(ParGridFunction()),
        _electrolyte_concentration_gf(ec),
        _electrolyte_concentration_gfc(&_electrolyte_concentration_gf),
        _scn(scn),
        _scp(scp),
        _jex_ne_tc(&_electrolyte_concentration_gfc, [=](real_t ec) { return kn * sqrt( _scn * ec * (1 - _scn) ); }),
        _jex_pe_tc(&_electrolyte_concentration_gfc, [=](real_t ec) { return kp * sqrt( _scp * ec * (1 - _scp) ); }) {}

      /// P2D
      ExchangeCurrentCoefficient(
        const real_t & kn,
        const real_t & kp,
        const ParGridFunction & sc,
        const ParGridFunction & ec):
        _surface_concentration_gf(sc),
        _electrolyte_concentration_gf(ec),
        _surface_concentration_gfc(&_surface_concentration_gf),
        _electrolyte_concentration_gfc(&_electrolyte_concentration_gf),
        _scn(0),
        _scp(0),
        _jex_ne_tc(&_surface_concentration_gfc, &_electrolyte_concentration_gfc, [=](real_t sc, real_t ec) { return kn * sqrt( sc * ec * (1 - sc) ); }),
        _jex_pe_tc(&_surface_concentration_gfc, &_electrolyte_concentration_gfc, [=](real_t sc, real_t ec) { return kp * sqrt( sc * ec * (1 - sc) ); }),
        _jex_pwc(Array<int>({NE, PE}), Array<Coefficient*>({static_cast<Coefficient*>(&_jex_ne_tc), static_cast<Coefficient*>(&_jex_pe_tc)})) {}

      /// SPM(e)
      virtual PWConstCoefficient & Eval()
      {
        /// SPMe
        if (_electrolyte_concentration_gfc.GetGridFunction())
        {
          ParFiniteElementSpace * x_fespace = _electrolyte_concentration_gf.ParFESpace();
          QuadratureSpace x_qspace(x_fespace->GetParMesh(), 2 * x_fespace->FEColl()->GetOrder());

          /// NE
          _jex_pwc.UpdateCoefficient(NE, _jex_ne_tc);
          real_t integral_ne = x_qspace.Integrate(_jex_pwc);
          _jex_pwc.ZeroCoefficient(NE);

          /// PE
          _jex_pwc.UpdateCoefficient(PE, _jex_pe_tc);
          real_t integral_pe = x_qspace.Integrate(_jex_pwc);
          _jex_pwc.ZeroCoefficient(PE);

          Vector c({integral_ne / NNE * NX, 0., integral_pe / NPE * NX});
          _jex_pwcc.UpdateConstants(c);
        }
        /// SPM
        else
        {
          _jex_pwcc(NE) = _jex_vec(NE - 1) * sqrt( _scn * (1 - _scn) );
          _jex_pwcc(PE) = _jex_vec(PE - 1) * sqrt( _scp * (1 - _scp) );
        }

        return _jex_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _jex_pwc.Eval(T, ip);
      }

};
