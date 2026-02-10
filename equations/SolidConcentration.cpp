#include "equations/SolidConcentration.hpp"

void SolidConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   MFEM_ASSERT(particle_region == NE || particle_region == PE, "Particle not in electrode!");

   // use of ternary operator to check if condition particle_region == NE is true
   // if true R is assigned to RN, else R is assigned to RP
   const real_t R = particle_region == NE ? RN : RP;
   const real_t D = particle_region == NE ? DN : DP;
   const real_t t_scale = particle_region == NE ? tn_scale : tp_scale;

   // Function Coefficient takes a lambda double and passes the physical coordinate of the evaluation point in vector &r
   // r(0) refers to the x coordinate and therefore the lambda returns x^2
   FunctionCoefficient r2([](const Vector & r){ return r(0) * r(0); });
   // pointwise product of two scalar coefficient (D/R^2 * r2)
   ProductCoefficient dr2(D / R / R, r2);
   // const_cast removes const from a references to a coeffcient
   ProductCoefficient jjr2(const_cast<Coefficient&>(j), r2);
   ProductCoefficient jr2(-1. / R / t_scale, jjr2);

   //delete the previous ParBilinearForm M and create a new bilinear form defined on the same finite element space
   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(r2));
   // Assemble(0) assembles the element matrices but does not finalize the global matrix
   M->Assemble(0); // keep sparsity pattern of M and K the same
   // applies the essential boundary conditions, finalize the global assembled matrix and produces a hypreParMatrix stored in Mmat
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   // loop over each element, evaluates dr2 at quadrature pt, compute shape function gradients and build the element stiffness matrix, then add it global sparse matrix
   K->AddDomainIntegrator(new DiffusionIntegrator(dr2));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   // surface_bdr defines which boundary attributes receive this integral
   Q->AddBoundaryIntegrator(new BoundaryLFIntegrator(jr2), const_cast<mfem::Array<int>&>(surface_bdr));
   Q->Assemble();
   // assemble the linear form into a parallel vector Qvec
   Qvec = std::move(*(Q->ParallelAssemble()));
   // zero out the essential dof for the right hand side vector
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   // b = K*x
   Kmat.Mult(x.GetBlock(SC + particle_id), b);
   // b = -K*x
   b.Neg();
   // b = -K*x + Q
   b += Qvec;
}

real_t SolidConcentration::SurfaceConcentration(const BlockVector &x)
{
   //if the current MPI rank owns the surface DOF, it will output a single DOF value, else it will return 0
   real_t csurf = IsSurfaceOwned() ? x.GetBlock(SC + particle_id)[surface_tdof] : 0;

   if (GetParticleRank() == GetSurfaceRank())
      return csurf;

   // communication send from particle mesh to surface mesh
   MPI_Request request;
   if (IsSurfaceOwned())
      MPI_Isend(&csurf, 1, MFEM_MPI_REAL_T, particle_rank, particle_id, MPI_COMM_WORLD, &request);
   if (IsParticleOwned())
      MPI_Recv(&csurf, 1, MFEM_MPI_REAL_T, surface_rank, particle_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   //ensures that send and recv is complete before carrying out next operation
   if (IsSurfaceOwned())
      MPI_Wait(&request, MPI_STATUS_IGNORE);

   return csurf;
}

int SolidConcentration::FindSurfaceTrueDof()
{
   Array<int> nat_dofs;
   // fills up nat_dofs with true DOF indices that lie on those boundaries
   fespace.GetEssentialTrueDofs(surface_bdr, nat_dofs);
   // This returns -1 if no DOFs are found on those boundaries, otherwise the first true DOF index on the surface
   return nat_dofs.IsEmpty() ? -1 : nat_dofs[0];
}

unsigned SolidConcentration::FindSurfaceRank()
{
   Array<bool> surface_rank(Mpi::WorldSize());
   // gather each rank's surface_owned status
   MPI_Allgather(&surface_owned, 1, MPI_CXX_BOOL, surface_rank.GetData(), 1, MPI_CXX_BOOL, MPI_COMM_WORLD);
   // returns the MPI rank ID where surface_owned == true
   return std::distance(surface_rank.begin(), std::find(surface_rank.begin(), surface_rank.end(), true));
}
void SolidConcentration::DebuggingCheck(const BlockVector &x)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(x.GetBlock(SC + particle_id));
   GridFunctionCoefficient u_gfc(&u_gf);
   FunctionCoefficient r2([](const Vector & x){ return x(0) * x(0); });
   ProductCoefficient ur2(u_gfc, r2);
   // quadrature order is taken to be 2 degree higher than the FE order for increased accuracy
   QuadratureSpace x_qspace(fespace.GetParMesh(), fespace.FEColl()->GetOrder() + 2);
   real_t integral = x_qspace.Integrate(ur2);

   if (!Mpi::WorldRank())
      std::cout << "Total flux accumulated (" << particle_id << ") = " << integral << std::endl;
}
