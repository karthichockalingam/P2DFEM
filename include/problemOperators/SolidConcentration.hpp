#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class SolidConcentration : public Equation
{
   private:
      size_t particle_id;

   public:
      SolidConcentration(ParFiniteElementSpace &f, size_t id) : Equation(f), particle_id(id) { nbc_bdr = 0; nbc_bdr[1] = 1; };
      virtual void update(const BlockVector &u, Coefficient &j = cc);
};
