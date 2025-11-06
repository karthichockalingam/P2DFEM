

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;


#pragma once

class SolidPotential : public Equation
{
   public:
      using Equation::Equation;
      virtual void Update(const BlockVector &x, const Coefficient &j, const real_t &dt);
};
