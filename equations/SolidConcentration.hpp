#include "mfem.hpp"
#include "Equation.hpp"

using namespace mfem;


#pragma once

class SolidConcentration : public Equation
{
   private:
      const unsigned particle_id;
      const unsigned particle_rank;
      const bool particle_owned;
      const int particle_dof;
      const Region particle_region;
      const Array<int> surface_bdr;
      const int surface_dof;
      const bool surface_owned;
      const unsigned surface_rank;

   public:
      SolidConcentration(ParFiniteElementSpace &f, const unsigned &id, const unsigned &rank, const int &dof, const Region &region)
         : Equation(f),
           particle_id(id), particle_rank(rank),
           particle_owned(particle_rank == Mpi::WorldRank()),
           particle_dof(dof), particle_region(region),
           surface_bdr({0, 1}),
           surface_dof(FindSurfaceDof()),
           surface_owned(surface_dof != -1),
           surface_rank(FindSurfaceRank())
      { }

      virtual void Update(const BlockVector &x, const Coefficient &j, const real_t &dt = 0.0);
      virtual real_t SurfaceConcentration(const BlockVector &x);
      bool IsParticleOwned(){ return particle_owned; }
      bool IsSurfaceOwned() { return surface_owned; }
      Region GetParticleRegion(){ return particle_region; }
      int GetParticleDof(){ return particle_dof; }
      unsigned GetSurfaceRank(){ return surface_rank; }
      int FindSurfaceDof();
      unsigned FindSurfaceRank();
};
