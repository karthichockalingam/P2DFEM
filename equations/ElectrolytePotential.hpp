

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      using Equation::Equation;
      virtual void update(const BlockVector &u, Coefficient &j = cc);
};
