#include "SolidConcentration.hpp"

void SolidConcentration::Update(const BlockVector &x, Coefficient &j)
{
   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   const real_t D = particle_region == PE ? DP : DN;
   const real_t A = particle_region == PE ? AP : AN;
   const real_t L = particle_region == PE ? LPE : LNE;
   const real_t S = particle_region == PE ? -1 : 1;

   FunctionCoefficient r2([](const Vector & r){ return r(0) * r(0); });
   FunctionCoefficient dr2([&](const Vector & r){ return D * r(0) * r(0); });
   FunctionCoefficient cr2([&](const Vector & r){ return S * I / A / L * r(0) * r(0); });

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(r2, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(dr2, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(cr2), nbc_bdr);
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

   real_t csurf = numeric_limits<real_t>::quiet_NaN();
   MPI_Request request;
   if (surface_ldof != -1)
      MPI_Isend(&r_gf[surface_ldof], 1, MFEM_MPI_REAL_T, particle_rank, 1, MPI_COMM_WORLD, &request);
   if (particle_rank == Mpi::WorldRank())
      MPI_Recv(&csurf, 1, MFEM_MPI_REAL_T, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   if (surface_ldof != -1)
      MPI_Wait(&request, MPI_STATUS_IGNORE);

   return csurf;
}
