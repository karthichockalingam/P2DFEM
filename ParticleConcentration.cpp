#include "ParticleConcentration.hpp"

const real_t D = 1.0;

double  function1(const Vector & x){ return D * x(0) * x(0); }
double  function2(const Vector & x){ return x(0) * x(0); }

void ParticleConcentration::update(const Vector &u)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   FunctionCoefficient diff(function1);
   FunctionCoefficient coeff(function2);
   ConstantCoefficient nbcCoef(1.0);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(coeff, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(diff, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(nbcCoef), nbc_bdr);
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   delete C;
   C = NULL; // re-compute C on the next ImplicitSolve
}
