#include "mfem.hpp"
using namespace std;
using namespace mfem;

#pragma once

class Equation
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
   Equation(ParFiniteElementSpace &f) : fespace(f), nbc_bdr(2), z(f.TrueVSize()) {};

   const HypreParMatrix & getM() const { return Mmat; };
   const HypreParMatrix & getK() const { return Kmat; };
   const HypreParVector & getQ() const { return Qvec; };

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   virtual void update(const BlockVector &u) = 0;

   virtual ~Equation()
   {
      delete C;
      delete M;
      delete K;
      delete Q;
   }
};