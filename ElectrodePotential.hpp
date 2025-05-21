

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrodePotential : public Equation
{
   public:
      using Equation::Equation;
      virtual void update(const BlockVector &u, GridFunction &ce, Coefficient &j = cc);
};


