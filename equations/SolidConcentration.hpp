#include "mfem.hpp"
#include "equations/Equation.hpp"

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
      const int surface_tdof;
      const bool surface_owned;
      const unsigned surface_rank;

   public:
      SolidConcentration(ParFiniteElementSpace &f, const unsigned &id, const unsigned &rank, const int &dof, const Region &region)
         : Equation(f),
           particle_id(id), particle_rank(rank),
           particle_owned(particle_rank == Mpi::WorldRank()),
           particle_dof(dof), particle_region(region),
           surface_bdr({0, 1}),
           surface_tdof(FindSurfaceTrueDof()),
           surface_owned(surface_tdof != -1),
           surface_rank(FindSurfaceRank())
      { }

      virtual void Update(const BlockVector &x, const Coefficient &j, const real_t &dt);
      virtual void Update(const BlockVector &, const GridFunctionCoefficient &, const Coefficient &, const real_t &) {}
      real_t SurfaceConcentration(const BlockVector &x);
      unsigned GetParticleRank(){ return particle_rank; }
      unsigned GetSurfaceRank(){ return surface_rank; }
      bool IsParticleOwned(){ return particle_owned; }
      bool IsSurfaceOwned() { return surface_owned; }
      Region GetParticleRegion(){ return particle_region; }
      int GetParticleDof(){ return particle_dof; }
      int FindSurfaceTrueDof();
      unsigned FindSurfaceRank();
      void DebuggingCheck(const BlockVector &x);
};
