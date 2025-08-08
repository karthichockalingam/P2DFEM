#include "SolidConcentration.hpp"

void SolidConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   const real_t R = particle_region == PE ? RP : RN;
   const real_t D = particle_region == PE ? DP : DN;
   const real_t t_scale = particle_region == PE ? tp_scale : tn_scale;

   FunctionCoefficient r2([](const Vector & r){ return r(0) * r(0); });
   ProductCoefficient dr2(D / R / R, r2);
   ProductCoefficient jjr2(const_cast<Coefficient&>(j), r2);
   ProductCoefficient jr2(-1. / R / t_scale, jjr2);

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
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(jr2), const_cast<mfem::Array<int>&>(surface_bdr));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(x.GetBlock(SC + particle_id), b);
   b.Neg();
   b += Qvec;
}

real_t SolidConcentration::SurfaceConcentration(const BlockVector &x)
{
   ParGridFunction r_gf(&fespace);
   r_gf.SetFromTrueDofs(x.GetBlock(SC + particle_id));

   real_t csurf = numeric_limits<real_t>::quiet_NaN();
   MPI_Request request;
   if (IsSurfaceOwned())
      MPI_Isend(&r_gf[surface_dof], 1, MFEM_MPI_REAL_T, particle_rank, 1, MPI_COMM_WORLD, &request);
   if (IsParticleOwned())
      MPI_Recv(&csurf, 1, MFEM_MPI_REAL_T, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   if (IsSurfaceOwned())
      MPI_Wait(&request, MPI_STATUS_IGNORE);

   return csurf;
}

int SolidConcentration::FindSurfaceDof()
{
   Array<int> nat_dofs;
   fespace.GetEssentialVDofs(surface_bdr, nat_dofs);
   return nat_dofs.Find(-1);
}