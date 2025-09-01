#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"
#include "coefficients/OpenCircuitPotentialCoefficient.hpp"
#include "coefficients/ReactionCurrentCoefficient.hpp"

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

   /// Gridfunctions defined over x_fespace (3 macro eqs plus _surface_ concentration)
   ParGridFunction * _ep_gf;
   ParGridFunction * _ec_gf;
   ParGridFunction * _sp_gf;
   ParGridFunction * _sc_gf;

   /// For the four gridfunctions over x_fespace (3 macro eqs plus _surface_ concentration)
   Array<int> block_offsets;
   /// For solution true vector (3 macros eqs plus NPAR _radial_ concentrations)
   Array<int> block_trueOffsets;
   /// For rhs true vectors (2 macro eqs)
   Array<int> potential_trueOffsets;
   /// For rhs true vectors (1 macro eq plus NPAR _radial_ concentrations)
   Array<int> concentration_trueOffsets;

   /// System matrices for concentration and potential eqs
   BlockOperator * Ac, * Ap;

   /// Block vector for the dofs of quantities defined over x_fespace (3 macro eqs plus _surface_ concentration)
   BlockVector _l;

   /// Reference to solution true dof vector
   BlockVector & _x;

   /// Constant timestep (as requested by the user, no time adaptivity supported)
   const real_t & _dt;

   /// Implicit solver for T = M + dt K
   CGSolver Solver;
   /// Preconditioner for the implicit solver
   HypreSmoother Prec;

   /// Auxiliary rhs vectors for concentrations and potential eqs
   mutable BlockVector bc, bp;

   /// File to write temporary data to
   std::ofstream file;

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
               const unsigned &ndofs, BlockVector &x, const real_t & dt);

   virtual void Mult(const Vector &x, Vector &dx_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   void SetGridFunctionsFromTrueVectors();

   virtual void Update();

   real_t GetSurfaceConcentration(const Region &r);
   void SetSurfaceConcentration();

   real_t ComputeReactionCurrent(const Region &r);
   ReactionCurrentCoefficient ComputeReactionCurrent();

   real_t ComputeExchangeCurrent(const Region &r);
   ExchangeCurrentCoefficient ComputeExchangeCurrent();

   real_t ComputeOpenCircuitPotential(const Region &r);
   OpenCircuitPotentialCoefficient ComputeOpenCircuitPotential();

   void ComputeVoltage(const BlockVector &x, real_t t, real_t dt);

   virtual void GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets);

   virtual ~P2DOperator()
   {
      delete ep;
      delete ec;
      delete sp;
      for (unsigned p = 0; p < NPAR; p++)
         delete sc[p];
      delete _ep_gf;
      delete _ec_gf;
      delete _sp_gf;
      delete _sc_gf;
      delete Ac;
      delete Ap;
      file.close();
   }
};
