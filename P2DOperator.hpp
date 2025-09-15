#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"
#include "coefficients/OpenCircuitPotentialCoefficient.hpp"
#include "coefficients/OverPotentialCoefficient.hpp"
#include "coefficients/ReactionCurrentCoefficient.hpp"

using namespace mfem;

#pragma once

class P2DOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace * &_x_fespace;
   Array<ParFiniteElementSpace *> &_r_fespace;

   ElectrolytePotential     *  _ep;
   SolidPotential           *  _sp;
   ElectrolyteConcentration *  _ec;
   Array<SolidConcentration *> _sc;

   /// Gridfunctions defined over x_fespace (3 macro eqs plus _surface_ concentration)
   ParGridFunction * _ep_gf;
   ParGridFunction * _sp_gf;
   ParGridFunction * _ec_gf;
   ParGridFunction * _sc_gf;

   /// For the four gridfunctions over _x_fespace (3 macro eqs plus _surface_ concentration)
   Array<int> _block_offsets;
   /// For solution true vector (3 macros eqs plus NPAR _radial_ concentrations)
   Array<int> _block_trueOffsets;
   /// For rhs true vectors (2 macro eqs)
   Array<int> _potential_trueOffsets;
   /// For rhs true vectors (1 macro eq plus NPAR _radial_ concentrations)
   Array<int> _concentration_trueOffsets;

   /// System matrices for concentration and potential eqs
   BlockOperator * _Ac, * _Ap;

   /// Block vector for the dofs of quantities defined over _x_fespace (3 macro eqs plus _surface_ concentration)
   BlockVector _l;

   /// Reference to solution true dof vector
   BlockVector & _x;

   /// Constant timestep (as requested by the user, no time adaptivity supported)
   const real_t & _dt;

   /// Implicit solver for T = M + dt K
   CGSolver _Solver;
   /// Preconditioner for the implicit solver
   HypreSmoother _Prec;

   /// Auxiliary rhs vectors for concentrations and potential eqs
   mutable BlockVector _bc, _bp;

   /// File to write temporary data to
   std::ofstream _file;

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
               const unsigned &ndofs, BlockVector &x, const real_t & dt);

   virtual void Mult(const Vector &x, Vector &dx_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   void SetGridFunctionsFromTrueVectors();

   virtual void UpdatePotentialEquations();
   virtual void UpdateConcentrationEquations();

   real_t GetSurfaceConcentration(const Region &r);
   void SetSurfaceConcentration();

   Array<real_t> ComputeParticleReactionCurrent();

   real_t ComputeReactionCurrent(const Region &r);
   ReactionCurrentCoefficient ComputeReactionCurrent();

   real_t ComputeExchangeCurrent(const Region &r);
   ExchangeCurrentCoefficient ComputeExchangeCurrent();

   real_t ComputeOpenCircuitPotential(const Region &r);
   OpenCircuitPotentialCoefficient ComputeOpenCircuitPotential();

   real_t ComputeOverPotential(const Region &r);
   OverPotentialCoefficient ComputeOverPotential();

   void ComputeVoltage(const BlockVector &x, real_t t, real_t dt);

   virtual void GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets);

   virtual ~P2DOperator()
   {
      delete _ep;
      delete _sp;
      delete _ec;
      for (unsigned p = 0; p < NPAR; p++)
         delete _sc[p];
      delete _ep_gf;
      delete _sp_gf;
      delete _ec_gf;
      delete _sc_gf;
      delete _Ac;
      delete _Ap;
      _file.close();
   }
};
