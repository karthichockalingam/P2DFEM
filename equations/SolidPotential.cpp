
#include "SolidPotential.hpp"
#include "../utils.hpp"

void SolidPotential::update(const BlockVector &u, Coefficient &j)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   real_t a = 1;

   ConstantCoefficient sigma(1.0);
   FunctionCoefficient dummy([](const Vector & x){ return sin(2*M_PI*x(0)); });
   ProductCoefficient source(dummy, j);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(sigma, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(u.GetBlock(SP), z);
   z.Neg();
   z += Qvec;
}
