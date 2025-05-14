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
   Solver.SetPreconditioner(Prec);

   const unsigned npar = r_fespace.Size();
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

   u.Update(block_trueOffsets);

   for (size_t p = 0; p < npar; p++)
      pc[p] = new ParticleConcentration(*r_fespace[p], p);
   
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //   M du_dt = -K(u + dt*du_dt) <=> (M + dt K) du_dt = -Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   
}
   
void P2DOperator::update(const Vector &u)
{
   // assemble B
   delete B;
   B = new BlockOperator(block_trueOffsets);
}
