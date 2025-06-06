#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class SolidConcentration : public Equation
{
   private:
      const size_t particle_id;
      const Region particle_region;
      const size_t particle_ltdof;

   public:
      SolidConcentration(ParFiniteElementSpace &f, size_t id, size_t ltdof = -1) : Equation(f), particle_id(id), particle_region(particle_id < NPEPAR ? PE : NE), particle_ltdof(ltdof) { nbc_bdr = 0; nbc_bdr[1] = 1; };
      virtual void update(const BlockVector &u, Coefficient &j = cc);
};
