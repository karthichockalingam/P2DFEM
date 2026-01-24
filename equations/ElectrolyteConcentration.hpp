

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;


#pragma once

class ElectrolyteConcentration : public Equation
{
   public:
      using Equation::Equation;
      virtual void Update(const BlockVector &, const Coefficient &, const real_t &) {}
      virtual void Update(const BlockVector &x, const GridFunctionCoefficient &ec_gfc, const Coefficient &j, const real_t &dt);
};
