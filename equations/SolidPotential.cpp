
#include "equations/SolidPotential.hpp"

void SolidPotential::Update(const BlockVector &x, const Coefficient &j)
{
   // Source term.
   // source_vec(x) = alpha * (L/N) * NX
   // source_vec = [S_NE, 0, S_PE]
   Vector source_vec({
      /* NE */ AN /* length scaling */ * (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ AP /* length scaling */ * (LPE / NPE * NX)});

   PWConstCoefficient source_part(source_vec);
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Effective conductivity (does not account for electrode filler).
   // sigma_vec(x) = (1-eps)*sigma*(L/N)*NX
   // sigma_vec = [sigma_NE, 0, sigma_PE]
   Vector sigma_vec({
      /* NE */ (1 - EPS_N) * SIGN /* length scaling */ / (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ (1 - EPS_P) * SIGP /* length scaling */ / (LPE / NPE * NX)});

   PWConstCoefficient sigma(sigma_vec);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(sigma));
   K->Assemble();
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   // assemble a parallel RHS vector with the essential boundary DOF zeroed out
   delete Q;
   Q = new ParLinearForm(&fespace);
   // Q = integral of (source(x) * phi(x))
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0);

   // b = Kx
   Kmat.Mult(x.GetBlock(SP), b);
   b.Neg();
   b += Qvec;
}
