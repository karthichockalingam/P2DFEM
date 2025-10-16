
#include "ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(x.GetBlock(EC));

   //real_t a = 1, tplus = 0;
   //ProductCoefficient source((1 - tplus) * a, const_cast<Coefficient&>(j));

   //ConstantCoefficient one(1.0);

   // Mass coefficient.
   Vector mass_vec({
      /* PE */  EPS_P /* length scaling */  * (LPE  / NPE  * NX),
      /* SEP */ EPS_S /* length scaling */  * (LSEP / NSEP * NX),
      /* NE */  EPS_N /* length scaling */  * (LNE  / NNE  * NX)});
   PWConstCoefficient mass_part(mass_vec);
   ConstantCoefficient t_scale(te_scale);
   ProductCoefficient mass(mass_part,t_scale);

   // Source term.
   Vector source_vec({
      /* PE */ (1 - TPLUS) * AP /* length scaling */ * (LPE / NPE * NX),
      /* SEP */ 0.,
      /* NE */ (1 - TPLUS) * AN /* length scaling */ * (LNE / NNE * NX)});

   PWConstCoefficient source_part(source_vec);
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Diffusion coefficient.
   Vector D_scale_vec({
      /* PE */  BPE  /* length scaling */ / pow(LPE / NPE * NX, 1) ,
      /* SEP */ BSEP  /* length scaling */ / pow(LSEP / NSEP * NX, 1)  ,
      /* NE */  BNE  /* length scaling */ / pow(LNE / NNE * NX, 1)  });

   DECoefficient D_coeff(u_gf);
   PWConstCoefficient D_scale_coeff(D_scale_vec);
   ProductCoefficient D(D_scale_coeff, D_coeff);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(mass));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(D));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(x.GetBlock(EC), b);
   b.Neg();
   b += Qvec;
}
