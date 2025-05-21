
#include "ElectrolytePotential.hpp"
#include "../utils.hpp"

void ElectrolytePotential::update(const BlockVector &u, Coefficient &j)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u.GetBlock(block_id - 2));

   real_t a = 1;

   ConstantCoefficient one(1.0);
   ConstantCoefficient kappa(1.0);
   ConstantCoefficient kappa_D(1.0);

   GridFunctionCoefficient EcCoeff(&u_gf);
   RatioCoefficient KappaOverEcCoeff(kappa_D, EcCoeff);
   GradientGridFunctionCoefficient GradEcCoeff(&u_gf);
   ScalarVectorProductCoefficient VecCoeff(KappaOverEcCoeff, GradEcCoeff);

   ProductCoefficient source(a, j);

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
   
   Kmat.Mult(u.GetBlock(block_id), z);
   z.Neg();
   z += Qvec;
}
