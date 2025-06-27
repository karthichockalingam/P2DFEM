#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"

using namespace std;
using namespace mfem;

#pragma once

class P2DOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace * &x_fespace;
   Array<ParFiniteElementSpace *> &r_fespace;

   ElectrolytePotential     * ep;
   ElectrolyteConcentration * ec;
   SolidPotential           * sp;
   Array<SolidConcentration *> sc;

   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   Array<int> nbc_bdr; // this list remains empty for pure Neumann b.c.

   Array<int> block_offsets;
   Array<int> block_trueOffsets;

   BlockOperator * A;

   HypreParMatrix *C; // C = M + dt K
   real_t current_dt;

   CGSolver Solver;    // Implicit solver for T = M + dt K
   HypreSmoother Prec; // Preconditioner for the implicit solver

   mutable BlockVector b; // auxiliary vector

   std::ofstream file; // file to write temporary data to

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
               const unsigned &ndofs, BlockVector &x);

   virtual void Mult(const Vector &x, Vector &dx_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   virtual void Update(const BlockVector &x);

   void ComputeVoltage(const BlockVector &x, real_t t, real_t dt);

   real_t ComputeExternalCurrent(const BlockVector &x);

   virtual void GetParticleLTDofs(Array<int> & particle_dofs, Array<int> & particle_offsets);

   virtual ~P2DOperator() { file.close(); }
};
