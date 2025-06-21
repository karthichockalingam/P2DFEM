#include "SolidConcentration.hpp"

void SolidConcentration::Update(const BlockVector &x, Coefficient &j)
{
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

   Kmat.Mult(x.GetBlock(SC + particle_id), b);
   b.Neg();
   b += Qvec;
}

real_t SolidConcentration::SurfaceConcentration(const BlockVector &x)
{
   Array<int> ess_ldofs;
   fespace.GetEssentialVDofs(nbc_bdr, ess_ldofs);

   assert(ess_ldofs.Sum() == -1 || ess_ldofs.Sum() == 0);

   int surface_ldof = ess_ldofs.Find(-1);

   ParGridFunction r_gf(&fespace);
   r_gf.SetFromTrueDofs(x.GetBlock(SC + particle_id));
   return (surface_ldof != -1) ? r_gf[surface_ldof] : numeric_limits<real_t>::quiet_NaN();
}
