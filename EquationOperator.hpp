#include "mfem.hpp"
using namespace std;
using namespace mfem;

#pragma once

class EquationOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   Array<int> nbc_bdr; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M;
   ParBilinearForm *K;
   ParLinearForm *Q;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParVector Qvec;

   HypreParMatrix *C; // C = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver C_solver;    // Implicit solver for T = M + dt K
   HypreSmoother C_prec; // Preconditioner for the implicit solver

   mutable Vector z; // auxiliary vector

public:
   EquationOperator(ParFiniteElementSpace &f, const Vector &u, const Array<int> &etl, const Array<int> &nb);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   virtual void SetParameters(const Vector &u) = 0;

   virtual ~EquationOperator() {}
};