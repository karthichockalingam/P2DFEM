#include "mfem.hpp"
#include "EquationOperator.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ParticleConcentration : public Equation
{
   public:
      using Equation::Equation;
      using Equation::~Equation;
      virtual void update(const Vector &u);
};