#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &u)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace), npar(r_fespace.Size()),
     B(NULL), current_dt(0.0), Solver(x_fespace->GetComm())
{
   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   //Solver.SetPreconditioner(Prec);

   const unsigned nb = 3 + npar; // 3 macro eqs + 1 micro eq/particle

   block_offsets.SetSize(nb + 1);
   block_trueOffsets.SetSize(nb + 1);
   block_offsets[0] = 0;
   block_trueOffsets[0] = 0;
   for (size_t p = 0; p < npar; p++)
   {
      block_offsets[p + 1] = r_fespace[p]->GetVSize();
      block_trueOffsets[p + 1] = r_fespace[p]->TrueVSize();
   }
   for (size_t x = 0; x < 3; x++)
   {
      block_offsets[npar + x + 1] = x_fespace->GetVSize();
      block_trueOffsets[npar + x + 1] = x_fespace->TrueVSize();
   }
   block_offsets.PartialSum();
   block_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << nb << std::endl;
      std::cout << "Unknowns (total): " << block_trueOffsets[nb] << std::endl;
   }

   u.Update(block_trueOffsets); u = 0.;
   z.Update(block_trueOffsets);

   for (size_t p = 0; p < npar; p++)
      pc.Append(new SolidConcentration(*r_fespace[p], p));
   ec = new ElectrolyteConcentration(*x_fespace, npar);
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //   M du_dt = -K(u + dt*du_dt) <=> (M + dt K) du_dt = -Ku
   // for du_dt, where K is linearized by using u from the previous timestep

   // assemble B
   for (size_t p = 0; p < npar; p++)
   {
      B->SetDiagonalBlock(p, Add(1, pc[p]->getM(), dt, pc[p]->getK()));
      z.GetBlock(p) = pc[p]->getZ();
   }

   B->SetDiagonalBlock(npar, Add(1, ec->getM(), dt, ec->getK()));
   z.GetBlock(npar) = ec->getZ();

   for (size_t x = 0; x < 2; x++)
   {
      IdentityOperator * I = new IdentityOperator(x_fespace->TrueVSize());
      B->SetDiagonalBlock(npar + x + 1, I);
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

   for (size_t p = 0; p < npar; p++)
      pc[p]->update(u);
   ec->update(u);
}
