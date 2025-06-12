#include "SolidConcentration.hpp"

const real_t D = 1.0;

real_t  function1(const Vector & x){ return D * x(0) * x(0); }
real_t  function2(const Vector & x){ return x(0) * x(0); }

void SolidConcentration::update(const BlockVector &u, Coefficient &j)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   const real_t D = particle_region == PE ? DP : DN;
   const real_t A = particle_region == PE ? AP : AN;
   const real_t L = particle_region == PE ? LPE : LNE;
   const real_t S = particle_region == PE ? -1 : 1;

   FunctionCoefficient x2([](const Vector & x){ return x(0) * x(0); });
   FunctionCoefficient dx2([&](const Vector & x){ return D * x(0) * x(0); });
   FunctionCoefficient cx2([&](const Vector & x){ return S * I / A / L * x(0) * x(0); });

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(x2, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(dx2, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(cx2), nbc_bdr);
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(u.GetBlock(SC + particle_id), z);
   z.Neg();
   z += Qvec;
}
