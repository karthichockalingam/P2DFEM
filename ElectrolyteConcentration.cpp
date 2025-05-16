
#include "ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::update(const BlockVector &u, Coefficient &j)
{
   ParGridFunction u_gf(&fespace);
   u_gf.SetFromTrueDofs(u);

   real_t a = 1, tplus = 0;
   ProductCoefficient source((1 - tplus) * a, j);

   IntegrationRule ir = IntRules.Get(Geometry::SEGMENT, 6);

   //std::cout << "Number of points " << ir.GetNPoints() << std::endl;

   ConstantCoefficient one(1.0);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(one, &ir));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(one, &ir));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?
   
   Kmat.Mult(u.GetBlock(block_id), z);
   z.Neg();
   z += Qvec;
}
