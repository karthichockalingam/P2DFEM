

#include "mfem.hpp"
#include "Equation.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolyteConcentration : public Equation
{
   public:
      using Equation::Equation;
      virtual void Update(const BlockVector &x, const Coefficient &j);
};
