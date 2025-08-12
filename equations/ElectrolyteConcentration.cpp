
#include "ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(x.GetBlock(EC));

   ConstantCoefficient one(1.0);

   // Mass coefficient.
   Vector mass_vec({
      /* PE */  EPS_P,               
      /* SEP */ EPS_S,
      /* NE */  EPS_N   });
   PWConstCoefficient mass(mass_vec);

   // Source term.
   Vector source_vec({
      /* PE */ (1 - TPLUS) * AP, 
      /* SEP */ 0.,
      /* NE */ (1 - TPLUS) * AN  });
   PWConstCoefficient source_part(source_vec);
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Diffusion coefficient.
   Vector D_scale_vec({
      /* PE */  De_p_scale, 
      /* SEP */ De_s_scale, 
      /* NE */  De_n_scale    });
   DECoefficient D_coeff(u_gf, ce_scale);
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
