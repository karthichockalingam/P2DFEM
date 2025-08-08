

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      using Equation::Equation;
      void Update(const BlockVector &x, const Coefficient &j, const real_t &dt);
};
