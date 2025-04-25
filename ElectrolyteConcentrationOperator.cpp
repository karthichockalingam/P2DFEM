
#include "ElectrolyteConcentrationOperator.hpp"

double  function3(const double & electrolyte_potential, const double & electrode_potential, const Vector & x)
{ 
   double val = (electrode_potential - electrolyte_potential) / (2.0);

   return std::sinh(val);
}

class FluxJGridFuncCoefficient : public Coefficient
 {
    GridFunction & _electrolyte_potential;
    GridFunction & _electrode_potential;
    function<double(const double &, const double &, const Vector &)>                GFunction;
    function<double(const double &, const double &, const Vector &, const double)>  TDGFunction;
 public:
 FluxJGridFuncCoefficient(GridFunction & electrolyte_potential, GridFunction & electrode_potential
                       , function<double(const double &, const double &, const Vector &)> foo)
                       : _electrolyte_potential(electrolyte_potential), _electrode_potential(electrode_potential), 
                       GFunction( move(foo) ) {};
 
                       FluxJGridFuncCoefficient(GridFunction & electrolyte_potential, GridFunction & electrode_potential
                       , function<double(const double &, const double &, const Vector &, const double)> foo)
                       : _electrolyte_potential(electrolyte_potential), _electrode_potential(electrode_potential),
                       TDGFunction( move(foo) ) {};
 
     double Eval(ElementTransformation &T, const IntegrationPoint &ip)
     {
       double x[3];
       Vector transip(x, 3);
       T.Transform(ip, transip);
       return (GFunction)? GFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), transip) 
                         : TDGFunction(_electrolyte_potential.GetValue(T, ip), _electrode_potential.GetValue(T, ip), transip, GetTime() );
     }
 };


void ElectrolyteConcentrationOperator::SetParameters(const Vector &u)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   real_t a, tplus;
   real_t temp = (1 - tplus) * a; 

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   ConstantCoefficient one(1.0);
   ConstantCoefficient nbcCoef(1.0);
   ConstantCoefficient Coeff(temp);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(one, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(one, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(Coeff));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   delete C;
   C = NULL; // re-compute C on the next ImplicitSolve
}

ElectrolyteConcentrationOperator::~ElectrolyteConcentrationOperator()
{
   delete C;
   delete M;
   delete K;
   delete Q;
}