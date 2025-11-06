
#include "equations/SolidPotential.hpp"

void SolidPotential::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   ConstantCoefficient sigma(1.0);
   FunctionCoefficient dummy([](const Vector & x){ return sin(2*M_PI*x(0)); });
   ProductCoefficient source(dummy, const_cast<Coefficient&>(j));

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(sigma));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(x.GetBlock(SP), b);
   b.Neg();
   b += Qvec;
   b *= 1./dt;
}
