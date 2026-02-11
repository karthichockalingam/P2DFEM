
#include "equations/ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::Update(const BlockVector &x, const GridFunctionCoefficient &ec_gfc, const Coefficient &j, const real_t &dt)
{
   // Mass coefficient.
   // mass_value  = epsilon * (L/N) * NX
   Vector mass_vec({
      /* NE */  EPS_N /* length scaling */  * (LNE  / NNE  * NX),
      /* SEP */ EPS_S /* length scaling */  * (LSEP / NSEP * NX),
      /* PE */  EPS_P /* length scaling */  * (LPE  / NPE  * NX)});
   // mass_part(x) = mass_vec[attribute(x)-1], attribute = 1,2,3   
   // creates a piecewise constant coefficient over the mesh
   PWConstCoefficient mass_part(mass_vec);
   ConstantCoefficient t_scale(te_scale);
   ProductCoefficient mass(mass_part, t_scale);

   // Source term = [value_NE, 0.0, value_PE].
   Vector source_vec({
      /* NE */ (1 - TPLUS) * AN /* length scaling */ * (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ (1 - TPLUS) * AP /* length scaling */ * (LPE / NPE * NX)});
   // creates a piecewise region dependent function
   PWConstCoefficient source_part(source_vec);
   //ProductCoefficient constructor takes a non-const coefficient, so constness must be removed
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Diffusion coefficient.
   // Diffusion_value = beta / (L/N) * NX
   Vector D_scale_vec({
      /* NE */  BNE  /* length scaling */ / (LNE  / NNE  * NX),
      /* SEP */ BSEP /* length scaling */ / (LSEP / NSEP * NX),
      /* PE */  BPE  /* length scaling */ / (LPE  / NPE  * NX)});

   // D_coeff(x) = DE(ex(x))
   TransformedCoefficient D_coeff(&const_cast<GridFunctionCoefficient&>(ec_gfc), [](real_t ec) { return DE(ec); });
   // D_scale(x) = D_scale_vec[attribute(x)-1], attribute=1,2,3
   PWConstCoefficient D_scale_coeff(D_scale_vec);
   // D(x) = D_scale(x) * DE(ec(x))
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
