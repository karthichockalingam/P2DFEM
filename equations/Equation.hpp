#include "mfem.hpp"
using namespace mfem;

#include "../constants.hpp"
using namespace constants;

using namespace std;

#pragma once

static ConstantCoefficient cc;

class Equation
{
protected:
   ParFiniteElementSpace &fespace;

   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   Array<int> nbc_bdr; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M = nullptr;
   ParBilinearForm *K = nullptr;
   ParLinearForm *Q = nullptr;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParVector Qvec;

   HypreSmoother prec; // Preconditioner for the implicit solver

   mutable Vector z; // auxiliary vector

public:
   Equation(ParFiniteElementSpace &f) : fespace(f), nbc_bdr(2), z(f.TrueVSize()) {};

   const HypreParMatrix & getM() const { return Mmat; };
   const HypreParMatrix & getK() const { return Kmat; };
   const Vector         & getZ() const { return z; };

   /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
   virtual void update(const BlockVector &u, Coefficient &j) = 0;

   virtual ~Equation()
   {
      delete M;
      delete K;
      delete Q;
   }
};
