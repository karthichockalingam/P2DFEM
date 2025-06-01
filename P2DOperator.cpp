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

   for (size_t p = 0; p < NPAR; p++)
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
   for (size_t p = 0; p < NPAR; p++)
      sc.Append(new SolidConcentration(*r_fespace[p], p, p < NPEPAR ? PE : NE));
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
   for (size_t p = 0; p < NPAR; p++)
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
   for (size_t p = 0; p < NPAR; p++)
      sc[p]->update(u);
}
