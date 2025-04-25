#include "mfem.hpp"
#include "EquationOperator.hpp"

using namespace std;
using namespace mfem;


#pragma once

class ParticleConcentrationOperator : public EquationOperator
{
   public:
      using EquationOperator::EquationOperator;
      virtual void SetParameters(const Vector &u);
      virtual ~ParticleConcentrationOperator();
};