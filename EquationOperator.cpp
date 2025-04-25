#include "EquationOperator.hpp"

EquationOperator::EquationOperator(ParFiniteElementSpace &f, const Vector &u, const Array<int> &etl, const Array<int> &nb)
   : TimeDependentOperator(f.GetTrueVSize(), (real_t) 0.0), fespace(f),
     ess_tdof_list(etl), nbc_bdr(nb), M(NULL), K(NULL), Q(NULL), C(NULL), current_dt(0.0),
     M_solver(f.GetComm()), C_solver(f.GetComm()), z(height)
{
   const real_t rel_tol = 1e-8;

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   C_solver.iterative_mode = false;
   C_solver.SetRelTol(rel_tol);
   C_solver.SetAbsTol(0.0);
   C_solver.SetMaxIter(100);
   C_solver.SetPrintLevel(0);
   C_solver.SetPreconditioner(C_prec);
}

void EquationOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   z += Qvec;
   HypreParMatrix A; Vector X, Z;
   ParGridFunction uu(&fespace), zz(&fespace);
   uu.SetFromTrueDofs(u);
   const SparseMatrix &R = *(fespace.GetRestrictionMatrix());
   R.MultTranspose(z, zz);
   K->FormLinearSystem(ess_tdof_list, uu, zz, A, X, Z);
   M_solver.Mult(Z, du_dt);
   du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void EquationOperator::ImplicitSolve(const real_t dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   if (!C)
   {
      C = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      C_solver.SetOperator(*C);
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
   C_solver.Mult(Z, du_dt);
   du_dt.SetSubVector(ess_tdof_list, 0.0);
}
