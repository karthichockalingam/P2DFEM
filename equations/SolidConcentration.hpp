#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class SolidConcentration : public Equation
{
   private:
      const unsigned particle_id;
      const unsigned particle_rank;
      const Region particle_region;
      const int particle_dof;

   public:
      SolidConcentration(ParFiniteElementSpace &f, unsigned id, unsigned rank, int dof = -1)
         : Equation(f),
           particle_id(id), particle_rank(rank),
           particle_region(id < NPEPAR ? PE : NE), particle_dof(dof)
      { nbc_bdr = 0; nbc_bdr[1] = 1; };
      virtual void Update(const BlockVector &x, Coefficient &j);
      virtual real_t SurfaceConcentration(const BlockVector &x);
      Region GetParticleRegion(){ return particle_region; }
      int GetParticleDof(){ return particle_dof; }
};
