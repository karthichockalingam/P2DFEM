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
      ParticleConcentration(ParFiniteElementSpace &f, size_t id) : Equation(f), id(id) {};
      virtual void update(const BlockVector &u);
};