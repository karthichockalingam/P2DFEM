
#include "ElectrolytePotential.hpp"
#include "../coefficients/utils.hpp"

void ElectrolytePotential::update(const BlockVector &u, Coefficient &j)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u.GetBlock(EC));

   real_t a = 1;

   ConstantCoefficient kappa(1.0);
   ConstantCoefficient kappa_D(1.0);
   FunctionCoefficient dummy([](const Vector & x){ return sin(2*M_PI*x(0)); });

   GridFunctionCoefficient EcCoeff(&u_gf);
   RatioCoefficient KappaOverEcCoeff(kappa_D, 1.0);
   GradientGridFunctionCoefficient GradEcCoeff(&u_gf);
   ScalarVectorProductCoefficient VecCoeff(KappaOverEcCoeff, GradEcCoeff);

   ProductCoefficient source(dummy, j);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(kappa, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->AddDomainIntegrator(new DomainLFGradIntegrator(VecCoeff));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   Kmat.Mult(u.GetBlock(EP), z);
   z.Neg();
   z += Qvec;
}
