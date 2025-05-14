#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ParticleConcentration : public Equation
{
   protected:
      size_t id;

   public:
      ParticleConcentration(ParFiniteElementSpace &f, size_t id) : Equation(f), id(id) { nbc_bdr = 0; nbc_bdr[1] = 1; };
      virtual void update(const BlockVector &u);
};
