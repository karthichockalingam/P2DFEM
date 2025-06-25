#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &x)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     A(NULL), current_dt(0.0), Solver(x_fespace->GetComm())
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

   x.Update(block_trueOffsets); x = 0.;
   b.Update(block_trueOffsets);

   ep = new ElectrolytePotential(*x_fespace);
   ec = new ElectrolyteConcentration(*x_fespace);
   sp = new SolidPotential(*x_fespace);

   if (M == SPM)
      for (unsigned p = 0; p < NPAR; p++)
         sc.Append(new SolidConcentration(*r_fespace[p], p, 0));
   else
   {
      Array<int> particle_dofs, particle_offsets;
      GetParticleLocalTrueDofs(particle_dofs, particle_offsets);
      for (unsigned p = 0; p < NPAR; p++)
      {
         unsigned rank = std::distance(particle_offsets.begin(),
            std::find_if(particle_offsets.begin(), particle_offsets.end(), [&](int i){ return p < i; }));
         unsigned offset = Mpi::WorldRank() > 0 ? particle_offsets[Mpi::WorldRank() - 1] : 0;
         bool mine = p >= offset && p < offset + particle_dofs.Size();
         int ltdof = mine ? particle_dofs[p - offset] : -1;
         sc.Append(new SolidConcentration(*r_fespace[p], p, rank, ltdof));
      }
   }
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &x, Vector &dx_dt)
{
   // Solve the equation:
   //   M dx_dt = -K(x + dt*dx_dt) <=> (M + dt K) dx_dt = -Kx
   // for dx_dt, where K is linearized by using x from the previous timestep

   // assemble A
   A->SetDiagonalBlock(EP, new HypreParMatrix(ep->GetK()));
   A->SetDiagonalBlock(EC, Add(1, ec->GetM(), dt, ec->GetK()));
   A->SetDiagonalBlock(SP, new HypreParMatrix(sp->GetK()));
   b.GetBlock(EP) = ep->GetZ();
   b.GetBlock(EC) = ec->GetZ();
   b.GetBlock(SP) = sp->GetZ();
   for (unsigned p = 0; p < NPAR; p++)
   {
      A->SetDiagonalBlock(SC + p, Add(1, sc[p]->GetM(), dt, sc[p]->GetK()));
      b.GetBlock(SC + p) = sc[p]->GetZ();
   }

   Solver.SetOperator(*A);
   Solver.Mult(b, dx_dt);
}

void P2DOperator::Update(const BlockVector &x)
{
   // rebuild A
   delete A;
   A = new BlockOperator(block_trueOffsets);
   A->owns_blocks = 1;

   // call point for j computation here

   ep->Update(x);
   ec->Update(x);
   sp->Update(x);
   for (unsigned p = 0; p < NPAR; p++)
      sc[p]->Update(x);

   for (unsigned p = 0; p < NPAR; p++)
   {
      real_t csurf = sc[p]->SurfaceConcentration(x);
      if (!isnan(csurf))
         std::cout << Mpi::WorldRank() << " " <<  p << " " << csurf << std::endl;
   }
}

void P2DOperator::GetParticleLocalTrueDofs(Array<int> & particle_dofs, Array<int> & particle_offsets)
{
   std::set<int> sep_global_dofs_set;
   for (int e = 0; e < x_fespace->GetNE(); e++)
      if (x_fespace->GetAttribute(e) == SEP)
      {
         Array<int> dofs;
         x_fespace->GetElementDofs(e, dofs);
         for (int d: dofs)
            sep_global_dofs_set.insert(x_fespace->GetGlobalTDofNumber(d));
      }

   unsigned max_sep_global_dofs = NSEP * x_fespace->GetElementOrder(0) + 1;
   Array<int> sep_global_dofs(max_sep_global_dofs); sep_global_dofs = -1;
   std::copy(sep_global_dofs_set.begin(), sep_global_dofs_set.end(), sep_global_dofs.begin());

   int all_sep_global_dofs[max_sep_global_dofs * Mpi::WorldSize()];
   MPI_Allgather(sep_global_dofs.GetData(), max_sep_global_dofs, MPI_INT,
                 all_sep_global_dofs, max_sep_global_dofs, MPI_INT, MPI_COMM_WORLD);

   sep_global_dofs_set = std::set<int>(all_sep_global_dofs, all_sep_global_dofs + max_sep_global_dofs * Mpi::WorldSize());
   
   std::set<int> particle_dofs_set;

   for (int d = 0; d < x_fespace->GetNDofs(); d++)
      if (sep_global_dofs_set.find(x_fespace->GetGlobalTDofNumber(d)) == sep_global_dofs_set.end())
         particle_dofs_set.insert(x_fespace->GetLocalTDofNumber(d));

   Array<int> boundary_dofs;
   x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (int d: boundary_dofs)
      particle_dofs_set.erase(d);
   particle_dofs_set.erase(-1);

   particle_dofs.SetSize(particle_dofs_set.size());
   std::copy(particle_dofs_set.begin(), particle_dofs_set.end(), particle_dofs.begin());

   int my_particles = particle_dofs.Size();
   particle_offsets.SetSize(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

   particle_offsets.PartialSum();
   assert(particle_offsets[Mpi::WorldSize() - 1] == NPAR);
}
