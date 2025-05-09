#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace &x_fespace, Array<ParFiniteElementSpace> &r_fespace,
                         const unsigned &ndofs, const Vector &u)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     B(NULL), current_dt(0.0), Solver(x_fespace.GetComm()), z(height)
{

   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   Solver.SetPreconditioner(Prec);
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   if (!C)
   {
      C = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      Solver.SetOperator(*B);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   z += Qvec;
   HypreParMatrix A; Vector X, Z;
   ParGridFunction uu(&fespace), zz(&fespace);
   uu.SetFromTrueDofs(u);
   const SparseMatrix &R = *(fespace.GetRestrictionMatrix());
   R.MultTranspose(z, zz);
   K->FormLinearSystem(ess_tdof_list, uu, zz, A, X, Z);
   Solver.Mult(Z, du_dt);
   du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void P2DOperator::update(const Vector &u)
{
   
}
