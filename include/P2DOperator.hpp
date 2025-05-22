#include "mfem.hpp"
#include "problemOperators/ElectrolytePotential.hpp"
#include "problemOperators/ElectrolyteConcentration.hpp"
#include "problemOperators/SolidPotential.hpp"
#include "problemOperators/SolidConcentration.hpp"

using namespace std;
using namespace mfem;

#pragma once

class P2DOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace * &x_fespace;
   Array<ParFiniteElementSpace *> &r_fespace;

   size_t npar;

   ElectrolytePotential     * ep;
   ElectrolyteConcentration * ec;
   SolidPotential           * sp;
   Array<SolidConcentration *> sc;

   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   Array<int> nbc_bdr; // this list remains empty for pure Neumann b.c.

   Array<int> block_offsets;
   Array<int> block_trueOffsets;

   BlockOperator * B;

   HypreParMatrix *C; // C = M + dt K
   real_t current_dt;

   CGSolver Solver;    // Implicit solver for T = M + dt K
   HypreSmoother Prec; // Preconditioner for the implicit solver

   mutable BlockVector z; // auxiliary vector

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
               const unsigned &ndofs, BlockVector &u);

   virtual void Mult(const Vector &u, Vector &du_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k) override;

   virtual void update(const BlockVector &u);

   virtual ~P2DOperator() {}
};
