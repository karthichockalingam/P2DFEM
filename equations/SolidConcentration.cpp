#include "SolidConcentration.hpp"

void SolidConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   const real_t R = particle_region == NE ? RN : RP;
   const real_t D = particle_region == NE ? DN : DP;
   const real_t t_scale = particle_region == NE ? tn_scale : tp_scale;

   FunctionCoefficient r2([](const Vector & r){ return r(0) * r(0); });
   ProductCoefficient dr2(D / R / R, r2);
   ProductCoefficient jjr2(const_cast<Coefficient&>(j), r2);
   ProductCoefficient jr2(-1. / R / t_scale, jjr2);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(r2));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(dr2));
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
   real_t csurf = IsSurfaceOwned() ? x.GetBlock(SC + particle_id)[surface_tdof] : 0;

   if (GetParticleRank() == GetSurfaceRank())
      return csurf;

   MPI_Request request;
   if (IsSurfaceOwned())
      MPI_Isend(&csurf, 1, MFEM_MPI_REAL_T, particle_rank, 1, MPI_COMM_WORLD, &request);
   if (IsParticleOwned())
      MPI_Recv(&csurf, 1, MFEM_MPI_REAL_T, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   if (IsSurfaceOwned())
      MPI_Wait(&request, MPI_STATUS_IGNORE);

   return csurf;
}

int SolidConcentration::FindSurfaceTrueDof()
{
   Array<int> nat_dofs;
   fespace.GetEssentialTrueDofs(surface_bdr, nat_dofs);
   return nat_dofs.IsEmpty() ? -1 : nat_dofs[0];
}

unsigned SolidConcentration::FindSurfaceRank()
{
   Array<bool> surface_rank(Mpi::WorldSize());
   MPI_Allgather(&surface_owned, 1, MPI_CXX_BOOL, surface_rank.GetData(), 1, MPI_CXX_BOOL, MPI_COMM_WORLD);
   return std::distance(surface_rank.begin(), std::find(surface_rank.begin(), surface_rank.end(), true));
}
void SolidConcentration::DebuggingCheck(const BlockVector &x)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(x.GetBlock(SC + particle_id));
   GridFunctionCoefficient u_gfc(&u_gf);
   FunctionCoefficient r2([](const Vector & x){ return x(0) * x(0); });
   ProductCoefficient ur2(u_gfc, r2);

   QuadratureSpace x_qspace(fespace.GetParMesh(), fespace.FEColl()->GetOrder() + 2);
   real_t integral = x_qspace.Integrate(ur2);

   if (!Mpi::WorldRank())
      std::cout << "Total flux accumulated (" << particle_id << ") = " << integral << std::endl;
}
