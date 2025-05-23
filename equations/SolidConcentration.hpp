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

   public:
      SolidConcentration(ParFiniteElementSpace &f, size_t id, Region r) : Equation(f), particle_id(id), particle_region(r) { nbc_bdr = 0; nbc_bdr[1] = 1; };
      virtual void update(const BlockVector &u, Coefficient &j = cc);
};
