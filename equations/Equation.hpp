#include "mfem.hpp"
using namespace mfem;

#include "../constants.hpp"
using namespace constants;

using namespace std;

#pragma once

class Equation
{
protected:
   ParFiniteElementSpace &fespace;

   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M = nullptr;
   ParBilinearForm *K = nullptr;
   ParLinearForm *Q = nullptr;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParVector Qvec;

   HypreSmoother prec; // Preconditioner for the implicit solver

   mutable Vector b; // auxiliary vector

public:
   Equation(ParFiniteElementSpace &f) : fespace(f), b(f.GetTrueVSize()) {};

   const HypreParMatrix & GetM() const { return Mmat; };
   const HypreParMatrix & GetK() const { return Kmat; };
   const Vector         & GetZ() const { return b; };

   /// Update the diffusion BilinearForm K using the given true-dof vector `x`.
   virtual void Update(const BlockVector &x, Coefficient &j) = 0;

   virtual ~Equation()
   {
      delete M;
      delete K;
      delete Q;
   }
};
