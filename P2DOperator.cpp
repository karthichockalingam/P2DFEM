#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &u)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     B(NULL), current_dt(0.0), Solver(x_fespace->GetComm())
{
   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   //Solver.SetPreconditioner(Prec);

   const unsigned nb = SC + NPAR; // 3 macro eqs + 1 micro eq/particle

   block_offsets.SetSize(nb + 1);
   block_trueOffsets.SetSize(nb + 1);

   block_offsets[0] = 0;
   block_offsets[EP + 1] = x_fespace->GetVSize();
   block_offsets[EC + 1] = x_fespace->GetVSize();
   block_offsets[SP + 1] = x_fespace->GetVSize();
   block_trueOffsets[0] = 0;
   block_trueOffsets[EP + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[EC + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[SP + 1] = x_fespace->GetTrueVSize();

   for (unsigned p = 0; p < NPAR; p++)
   {
      block_offsets[SC + p + 1] = r_fespace[p]->GetVSize();
      block_trueOffsets[SC + p + 1] = r_fespace[p]->GetTrueVSize();
   }

   block_offsets.PartialSum();
   block_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << nb << std::endl;
      std::cout << "Unknowns (rank 0): " << block_trueOffsets[nb] << std::endl;
   }

   u.Update(block_trueOffsets); u = 0.;
   z.Update(block_trueOffsets);

   ep = new ElectrolytePotential(*x_fespace);
   ec = new ElectrolyteConcentration(*x_fespace);
   sp = new SolidPotential(*x_fespace);

   if (M == SPM)
      for (unsigned p = 0; p < NPAR; p++)
         sc.Append(new SolidConcentration(*r_fespace[p], p));
   else
   {
      Array<int> particle_dofs; unsigned particle_offset;
      GetParticleLocalTrueDofs(particle_dofs, particle_offset);
      for (unsigned p = 0; p < NPAR; p++)
      {
         bool my_particle = p >= particle_offset && p < particle_offset + particle_dofs.Size();
         sc.Append(new SolidConcentration(*r_fespace[p], p, my_particle ? particle_dofs[p - particle_offset] : -1));
      }
   }
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //   M du_dt = -K(u + dt*du_dt) <=> (M + dt K) du_dt = -Ku
   // for du_dt, where K is linearized by using u from the previous timestep

   // assemble B
   B->SetDiagonalBlock(EP, new HypreParMatrix(ep->getK()));
   B->SetDiagonalBlock(EC, Add(1, ec->getM(), dt, ec->getK()));
   B->SetDiagonalBlock(SP, new HypreParMatrix(sp->getK()));
   z.GetBlock(EP) = ep->getZ();
   z.GetBlock(EC) = ec->getZ();
   z.GetBlock(SP) = sp->getZ();
   for (unsigned p = 0; p < NPAR; p++)
   {
      B->SetDiagonalBlock(SC + p, Add(1, sc[p]->getM(), dt, sc[p]->getK()));
      z.GetBlock(SC + p) = sc[p]->getZ();
   }

   Solver.SetOperator(*B);
   Solver.Mult(z, du_dt);
}

void P2DOperator::update(const BlockVector &u)
{
   // rebuild B
   delete B;
   B = new BlockOperator(block_trueOffsets);
   B->owns_blocks = 1;

   // call point for j computation here

   ep->update(u);
   ec->update(u);
   sp->update(u);
   for (unsigned p = 0; p < NPAR; p++)
      sc[p]->update(u);
}

void P2DOperator::GetParticleLocalTrueDofs(Array<int> & particle_dofs, unsigned & particle_offset)
{
   std::set<int> particle_dofs_set;

   for (int d = 0; d < x_fespace->GetNDofs(); d++)
      particle_dofs_set.insert(x_fespace->GetLocalTDofNumber(d));

   for (int e = 0; e < x_fespace->GetNE(); e++)
   {
      if (x_fespace->GetAttribute(e) != SEP)
         continue;

      Array<int> dofs;
      x_fespace->GetElementDofs(e, dofs);
      for (int d: dofs)
         particle_dofs_set.erase(x_fespace->GetLocalTDofNumber(d));
   }

   Array<int> boundary_dofs;
   x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (int d: boundary_dofs)
      particle_dofs_set.erase(d);
   particle_dofs_set.erase(-1);

   particle_dofs.SetSize(particle_dofs_set.size());
   std::copy(particle_dofs_set.begin(), particle_dofs_set.end(), particle_dofs.begin());

   HYPRE_BigInt my_particles = particle_dofs.Size();
   Array<HYPRE_BigInt> all_particles(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, HYPRE_MPI_BIG_INT,
                 all_particles.GetData(), 1, HYPRE_MPI_BIG_INT, MPI_COMM_WORLD);

   all_particles.PartialSum();
   particle_offset = Mpi::WorldRank() > 0 ? all_particles[Mpi::WorldRank() - 1] : 0;

   // FIX ME: Guarantees we crash in cases the partition boundary coincides with the electrode/separator boundary
   assert(all_particles[Mpi::WorldSize() - 1] == NPAR);
}
