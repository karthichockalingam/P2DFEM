

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      using Equation::Equation;
      virtual void Update(const BlockVector &x, Coefficient &j = cc);
};
