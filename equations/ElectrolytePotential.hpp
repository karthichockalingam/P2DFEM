

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      using Equation::Equation;
      void Update(const BlockVector &x, const Coefficient &j, const real_t &dt);
};
