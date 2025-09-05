#include "mfem.hpp"
using namespace mfem;

class ExchangeCurrentCoefficient: public Coefficient
{
   private:
      const ParGridFunction & _surface_concentration_gf;
      const ParGridFunction & _electrolyte_concentration_gf;

      GridFunctionCoefficient _surface_concentration_gfc;
      GridFunctionCoefficient _electrolyte_concentration_gfc;

      TransformedCoefficient _jex_ne_tc;
      ConstantCoefficient _jex_sep_cc;
      TransformedCoefficient _jex_pe_tc;

      Vector _jex_vec;
      PWConstCoefficient _jex_pwcc;
      PWCoefficient _jex_pwc;
      Coefficient & _jex;

   public:
      /// Default
      ExchangeCurrentCoefficient():
        _surface_concentration_gf(ParGridFunction()),
        _electrolyte_concentration_gf(ParGridFunction()),
        _jex_ne_tc(nullptr, nullptr, [](real_t, real_t) { return 0; }),
        _jex_pe_tc(nullptr, nullptr, [](real_t, real_t) { return 0; }),
        _jex(_jex_pwcc) {}

      /// SPM
      ExchangeCurrentCoefficient(
        const real_t & kn,
        const real_t & kp,
        const real_t & scn,
        const real_t & scp,
        const real_t & ec):
        _surface_concentration_gf(ParGridFunction()),
        _electrolyte_concentration_gf(ParGridFunction()),
        _jex_ne_tc(nullptr, nullptr, [](real_t, real_t) { return 0; }),
        _jex_pe_tc(nullptr, nullptr, [](real_t, real_t) { return 0; }),
        _jex_vec({kn * sqrt( scn * ec * (1 - scn) ), 0., kp * sqrt( scp * ec * (1 - scp) )}),
        _jex_pwcc(_jex_vec),
        _jex(_jex_pwcc) {}

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
        _jex_ne_tc(&_electrolyte_concentration_gfc, [=](real_t ec) { return kn * sqrt( scn * ec * (1 - scn) ); }),
        _jex_sep_cc(0),
        _jex_pe_tc(&_electrolyte_concentration_gfc, [=](real_t ec) { return kp * sqrt( scp * ec * (1 - scp) ); }),
        _jex(_jex_pwcc) {}

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
        _jex_ne_tc(&_surface_concentration_gfc, &_electrolyte_concentration_gfc, [=](real_t sc, real_t ec) { return kn * sqrt( sc * ec * (1 - sc) ); }),
        _jex_sep_cc(0),
        _jex_pe_tc(&_surface_concentration_gfc, &_electrolyte_concentration_gfc, [=](real_t sc, real_t ec) { return kp * sqrt( sc * ec * (1 - sc) ); }),
        _jex_pwc(Array<int>({NE, SEP, PE}), Array<Coefficient*>({static_cast<Coefficient*>(&_jex_ne_tc), static_cast<Coefficient*>(&_jex_sep_cc), static_cast<Coefficient*>(&_jex_pe_tc)})),
        _jex(_jex_pwc) {}

      /// SPM(e)
      virtual PWConstCoefficient Eval()
      {
        MFEM_ASSERT(&_jex == &_jex_pwcc, "ExchangeCurrentCoefficient does not wrap a PWConstCoefficient");

        /// SPMe
        if (!_jex_pwcc.GetNConst())
        {
          ParFiniteElementSpace * x_fespace = _electrolyte_concentration_gf.ParFESpace();
          Array<int> markers(x_fespace->GetParMesh()->attributes.Max());

          // NE
          markers = 0; markers[NE - 1] = 1;
          ParLinearForm sum_ne(x_fespace);
          sum_ne.AddDomainIntegrator(new DomainLFIntegrator(_jex_ne_tc), markers);
          sum_ne.Assemble();

          // PE
          markers = 0; markers[PE - 1] = 1;
          ParLinearForm sum_pe(x_fespace);
          sum_pe.AddDomainIntegrator(new DomainLFIntegrator(_jex_pe_tc), markers);
          sum_pe.Assemble();

          Vector c({x_fespace->GetParMesh()->ReduceInt(sum_ne.Sum()) / LNE,
                    0.,
                    x_fespace->GetParMesh()->ReduceInt(sum_pe.Sum()) / LPE});
          _jex_pwcc.UpdateConstants(c);
        }

        return _jex_pwcc;
      }

      /// P2D
      virtual real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
        return _jex.Eval(T, ip);
      }

};
