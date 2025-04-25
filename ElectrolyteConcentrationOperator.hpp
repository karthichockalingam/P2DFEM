

#include "mfem.hpp"
#include "EquationOperator.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ElectrolyteConcentrationOperator : public EquationOperator
{
   public:
      using EquationOperator::EquationOperator;
      virtual void SetParameters(const Vector &u);
      virtual ~ElectrolyteConcentrationOperator();
};