#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"

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

   Array<int> block_trueOffsets;
   Array<int> concentration_trueOffsets;
   Array<int> potential_trueOffsets;

   BlockOperator * Ac, * Ap; // system matrices (LHS)

   real_t current_dt;

   CGSolver Solver;    // Implicit solver for T = M + dt K
   HypreSmoother Prec; // Preconditioner for the implicit solver

   mutable BlockVector bc, bp; // auxiliary vector (RHS)

   std::ofstream file; // file to write temporary data to

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
               const unsigned &ndofs, BlockVector &x);

   virtual void Mult(const Vector &x, Vector &dx_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   virtual void Update(const BlockVector &x, const real_t &dt);

   PWConstCoefficient ComputeReactionCurrent();
   ConstantCoefficient ComputeReactionCurrent(const Region &r);
   ConstantCoefficient ComputeReactionCurrent(const BlockVector &x);

   real_t GetSurfaceConcentration(const Region &r, const BlockVector &x);
   ParGridFunction GetSurfaceConcentration(const BlockVector &x);

   real_t ComputeExchangeCurrent(const Region &r, const BlockVector &x);
   ExchangeCurrentCoefficient ComputeExchangeCurrent(const BlockVector &x);

   real_t ComputeOpenCircuitPotential(const Region &r, const real_t &x);

   void ComputeVoltage(const BlockVector &x, real_t t, real_t dt);

   virtual void GetParticleDofs(Array<int> & particle_dofs, Array<int> & particle_offsets);

   virtual ~P2DOperator() { file.close(); }
};
