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
      const bool particle_owned;
      const Region particle_region;
      const int particle_dof;
      const Array<int> surface_bdr;
      const int surface_dof_rank;
      const int surface_dof;
      const bool surface_owned;

   public:
      SolidConcentration(ParFiniteElementSpace &f, unsigned id, unsigned rank, unsigned surface_rank, int dof = -1)
         : Equation(f),
           particle_id(id), particle_rank(rank),
           particle_owned(particle_rank == Mpi::WorldRank()),
           particle_region(id < NPEPAR ? PE : NE), particle_dof(dof),
           surface_bdr({0, 1}), 
           surface_dof_rank(surface_rank),
           surface_dof(FindSurfaceDof()),
           surface_owned(surface_dof != -1)
      { }

      virtual void Update(const BlockVector &x, const Coefficient &j, const real_t &dt = 0.0);
      virtual real_t SurfaceConcentration(const BlockVector &x);
      bool IsParticleOwned(){ return particle_owned; }
      bool IsSurfaceOwned() { return surface_owned; }
      Region GetParticleRegion(){ return particle_region; }
      int GetParticleDof(){ return particle_dof; }
      int GetSurfaceDofRank() {return surface_dof_rank;}
      int FindSurfaceDof();
};
