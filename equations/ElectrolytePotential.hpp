

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      using Equation::Equation;
      void Update(const BlockVector &x, const Coefficient &j, real_t dt);
};
