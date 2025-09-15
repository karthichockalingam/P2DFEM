
#include "ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(x.GetBlock(EC));

   ConstantCoefficient one(1.0);
   ConstantCoefficient neg_one(-1.0);

   // Mass coefficient.
   Vector mass_vec({
      /* PE */  EPS_P * /* length scaling */ (LPE  / NPE  * NX),
      /* SEP */ EPS_S * /* length scaling */ (LSEP / NSEP * NX),
      /* NE */  EPS_N * /* length scaling */ (LNE  / NNE  * NX)});
   PWConstCoefficient mass_part(mass_vec);
   ConstantCoefficient t_scale(te_scale);
   //ConstantCoefficient t_scale(1.0);
   ProductCoefficient mass(mass_part,t_scale);

   /*std::cout << " " << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << "MASS TERMS:" << std::endl;
   std::cout << "param.NE.eps = " << EPS_N << std::endl;
   std::cout << "param.SEP.eps = " << EPS_S << std::endl;
   std::cout << "param.PE.eps = " << EPS_P << std::endl;
   std::cout << "te_scale = " << te_scale << std::endl;
   std::cout << " " << std::endl;
   std::cout << "Mass (PE, SEP, NE) = " << mass_vec(0)*te_scale << ", " << mass_vec(1)*te_scale << ", " << mass_vec(2)*te_scale << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << " " << std::endl;*/

   // Source term.
   Vector source_vec({
      /* PE */ (1 - TPLUS) * AP * /* length scaling */ (LPE / NPE * NX),
      /* SEP */ 0.,
      /* NE */ (1 - TPLUS) * AN * /* length scaling */ (LNE / NNE * NX)});

   Vector j_vec({
      /* PE */ -5.980666,
      /* SEP */ 0.,
      /* NE */ 5.2822535  });

   PWConstCoefficient source_part(source_vec);
   PWConstCoefficient j_vec_coeff(j_vec);
   //ProductCoefficient source(source_part, const_cast<Coefficient&>(j));
   ProductCoefficient source(source_part, j_vec_coeff);
   //ProductCoefficient source_neg(source, neg_one);


   /*std::cout << " " << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << "SOURCE TERMS:" << std::endl;
   std::cout << "tplus = " << TPLUS << std::endl;
   std::cout << "param.NE.as = " << AN << std::endl;
   std::cout << "param.PE.as = " << AP << std::endl;
   std::cout << "(1 - TPLUS) * AN = " << (1 - TPLUS) * AN << std::endl;
   std::cout << "(1 - TPLUS) * AP = " << (1 - TPLUS) * AP << std::endl;
   //std::cout << "jn = " << j_n << std::endl;
   //std::cout << "jp = " << j_p << std::endl;
   std::cout << " " << std::endl;
   std::cout << "Source (PE, SEP, NE) = " << source_vec(0) * -5.980666 << ", " << source_vec(1) << ", " << source_vec(2) * 5.2822535 << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << " " << std::endl;*/

   // Diffusion coefficient.
   Vector D_scale_vec({
      /* PE */  De_p_scale / /* length scaling */ (LPE / NPE * NX),
      /* SEP */ De_s_scale / /* length scaling */ (LSEP / NSEP * NX),
      /* NE */  De_n_scale / /* length scaling */ (LNE / NNE * NX)});

   DECoefficient D_coeff(u_gf, ce_scale);
   //ConstantCoefficient D_coeff_fixed(0.5);
   FunctionCoefficient D_coeff_func([](const Vector & x){ return DE(1.0); });
   PWConstCoefficient D_scale_coeff(D_scale_vec);
   //ProductCoefficient D(D_scale_coeff, D_coeff_func);
   ProductCoefficient D(D_scale_coeff, D_coeff);
   //ProductCoefficient negD(D,neg_one);

   //std::cout << "D_scale_vec = " << D_scale_vec(0) << ", " << D_scale_vec(1) << ", " << D_scale_vec(2) << std::endl;
   //std::cout << "D(1.5 * ce_scale) = " << DE(1.5 * ce_scale) << std::endl;
   /*std::cout << " " << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << "DIFFUSION TERMS:" << std::endl;
   std::cout << "D(1. * ce_scale) * De_p_scale = " << DE(1.) * De_p_scale << std::endl;
   std::cout << "D(1. * ce_scale) * De_s_scale = " << DE(1.) * De_s_scale << std::endl;
   std::cout << "D(1. * ce_scale) * De_n_scale = " << DE(1.) * De_n_scale << std::endl;
   std::cout << "****************" << std::endl;
   std::cout << " " << std::endl;*/

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

   Mmat.Print("Mmat.txt");
   Kmat.Print("Kmat.txt");

   //mfem_error("STOP");

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(x.GetBlock(EC), b);
   b.Neg();
   b += Qvec;

   /*std::cout << "b:" << endl;
   for (int i = 0; i < b.Size(); i++)
   {
      std::cout << b[i] << std::endl;
   }

   std::cout << "Qvec:" << std::endl;
   for (int i = 0; i < Qvec.Size(); i++)
   {
      std::cout << Qvec[i] << std::endl;
   }*/
}
