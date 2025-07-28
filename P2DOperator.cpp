#include "P2DOperator.hpp"
#include "utils.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &x)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     A(NULL), current_dt(0.0), Solver(x_fespace->GetComm()), file("data.csv")
{
   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   //Solver.SetPreconditioner(Prec);

   const unsigned nb = SC + NPAR; // 3 macro eqs + 1 micro eq/particle

   block_offsets.SetSize(nb + 1);
   block_trueOffsets.SetSize(nb + 1);

   block_offsets[0] = 0;
   block_offsets[EP + 1] = x_fespace->GetVSize();
   block_offsets[EC + 1] = x_fespace->GetVSize();
   block_offsets[SP + 1] = x_fespace->GetVSize();
   block_trueOffsets[0] = 0;
   block_trueOffsets[EP + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[EC + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[SP + 1] = x_fespace->GetTrueVSize();

   for (unsigned p = 0; p < NPAR; p++)
   {
      block_offsets[SC + p + 1] = r_fespace[p]->GetVSize();
      block_trueOffsets[SC + p + 1] = r_fespace[p]->GetTrueVSize();
   }

   block_offsets.PartialSum();
   block_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << nb << std::endl;
      std::cout << "Unknowns (rank 0): " << block_trueOffsets[nb] << std::endl;
   }

   x.Update(block_trueOffsets); x = 0.;
   b.Update(block_trueOffsets);

   ep = new ElectrolytePotential(*x_fespace);
   ec = new ElectrolyteConcentration(*x_fespace);
   sp = new SolidPotential(*x_fespace);

   if (M == SPM)
      for (unsigned p = 0; p < NPAR; p++)
         sc.Append(new SolidConcentration(*r_fespace[p], p, 0));
   else
   {
      Array<int> particle_dofs, particle_offsets;
      GetParticleDofs(particle_dofs, particle_offsets);
      for (unsigned p = 0; p < NPAR; p++)
      {
         unsigned rank = std::distance(particle_offsets.begin(),
            std::find_if(particle_offsets.begin(), particle_offsets.end(), [&](int i){ return p < i; }));
         unsigned offset = Mpi::WorldRank() > 0 ? particle_offsets[Mpi::WorldRank() - 1] : 0;
         bool mine = p >= offset && p < offset + particle_dofs.Size();
         int dof = mine ? particle_dofs[p - offset] : -1;
         sc.Append(new SolidConcentration(*r_fespace[p], p, rank, dof));
      }
   }
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &x, Vector &dx_dt)
{
   // Solve the equation:
   //   M dx_dt = -K(x + dt*dx_dt) <=> (M + dt K) dx_dt = -Kx
   // for dx_dt, where K is linearized by using x from the previous timestep

   // assemble A
   A->SetDiagonalBlock(EP, new HypreParMatrix(ep->GetK()));
   A->SetDiagonalBlock(EC, Add(1, ec->GetM(), dt, ec->GetK()));
   A->SetDiagonalBlock(SP, new HypreParMatrix(sp->GetK()));
   b.GetBlock(EP) = ep->GetZ();
   b.GetBlock(EC) = ec->GetZ();
   b.GetBlock(SP) = sp->GetZ();
   for (unsigned p = 0; p < NPAR; p++)
   {
      A->SetDiagonalBlock(SC + p, Add(1, sc[p]->GetM(), dt, sc[p]->GetK()));
      b.GetBlock(SC + p) = sc[p]->GetZ();
   }

   Solver.SetOperator(*A);
   Solver.Mult(b, dx_dt);
}

void P2DOperator::Update(const BlockVector &x)
{
   // rebuild A
   delete A;
   A = new BlockOperator(block_trueOffsets);
   A->owns_blocks = 1;

   FunctionCoefficient j = ComputeReactionCurrent(x);

   ep->Update(x, j);
   ec->Update(x, j);
   sp->Update(x, j);
   for (unsigned p = 0; p < NPAR; p++)
      sc[p]->Update(x, j);

}

FunctionCoefficient P2DOperator::ComputeReactionCurrent(const BlockVector &x)
{
   const real_t jp = - I / AP / LPE;
   const real_t jn = + I / AN / LNE;
   auto j = [=](const Vector & p){ return p(0) < LPE ? jp : p(0) < LPE + LSEP ? 0 : jn; };

   std::cout << "jp = " << jp << std::endl;
   std::cout << "jn = " << jn << std::endl;
   std::cout << "AP = " << AP << std::endl;
   std::cout << "AN = " << AN << std::endl;

   return FunctionCoefficient(j);
}

ConstantCoefficient P2DOperator::ComputeExchangeCurrent(const BlockVector &x)
{
   ParGridFunction cs_gf(x_fespace);
   cs_gf = 0;

   for (unsigned p = 0; p < NPAR; p++)
   {
      real_t csurf = sc[p]->SurfaceConcentration(x);
      if (!isnan(csurf))
         cs_gf(sc[p]->GetParticleDof()) = csurf;
   }
   // Apply prolongation after restriction. Might be unnecessary, but guarantees
   // all processors have the right information for all their local dofs. 
   cs_gf.SetFromTrueVector();

   ParGridFunction ec_gf(x_fespace);
   ec_gf.SetFromTrueDofs(x.GetBlock(EC));

   ExchangeCurrentCoefficient coeff(cs_gf, ec_gf);
   ParLinearForm sum(x_fespace);
   sum.AddDomainIntegrator(new DomainLFIntegrator(coeff));
   sum.Assemble();

   real_t reduction_result = sum.Sum();
   MPI_Allreduce(MPI_IN_PLACE, &reduction_result, 1, MFEM_MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);
   //Is it correct to use the total length or is it the lenght of the region?
   return ConstantCoefficient(reduction_result / (LPE + LSEP + LNE));
}

void P2DOperator::ComputeVoltage(const BlockVector &x, real_t t, real_t dt)
{

   real_t csurf[NPAR];
   for (unsigned p = 0; p < NPAR; p++)
   {
      csurf[p] = sc[p]->SurfaceConcentration(x);
      if (!isnan(csurf[p]))
         std::cout << "[Rank " << Mpi::WorldRank() << "]"
                   << " Surface concentration (" << p << ") = "
                   << csurf[p] << std::endl;
   }

   //real_t voltage = 10 - csurf[1]/10 +
   //             asinh(- I / AP / LPE / 2 / sqrt((10+csurf[0])*-csurf[0])) -
   //             asinh(  I / AN / LNE / 2 / sqrt(csurf[1]*(10-csurf[1])));

   real_t cp0 = 0.27;
   real_t cn0 = 0.94;

   real_t cp = csurf[0] + cp0;   // Particle surface concentration at the positive electrode.
   real_t cn = csurf[1] + cn0;   // Particle surface concentration at the negative electrode.

   // Definition from LIONSIMBA: https://doi.org/10.1149/2.0291607jes
   real_t theta_p = cp; // As cp is non-dimensionalised, theta_p = cp.
   real_t theta_n = cn; // As cn is non-dimensionalised, theta_n = cn.

   std::cout << "DP = " << DP << std::endl;
   std::cout << "DN = " << DN << std::endl;

   // Open Circuit Potential (no temperature dependence).
   real_t Up = -0.8090*theta_p + 4.4875 - 0.0428*tanh(18.5138*(theta_p - 0.5542)) - 17.7326*tanh(15.7890*(theta_p - 0.3117)) + 17.5842*tanh(15.9308*(theta_p - 0.3120));
   real_t Un = 1.97938*exp(-39.3631*theta_n) + 0.2482 - 0.0909*tanh(29.8538*(theta_n - 0.1234)) - 0.04478*tanh(14.9159*(theta_n - 0.2769)) - 0.0205*tanh(30.4444*(theta_n - 0.6103));


   real_t T = 1.0;

   real_t jp = - I / AP / LPE;
   real_t jn = + I / AN / LNE;

   real_t ce0_scaled = 1.0;

   real_t j0_p =  KP * pow(cp * ce0_scaled * abs(1.0 - cp), 0.5);
   real_t j0_n =  KN * pow(cn * ce0_scaled * abs(1.0 - cn), 0.5);

   real_t eta_p = 2 * T * asinh(jp / 2.0 / j0_p);
   real_t eta_n = 2 * T * asinh(jn / 2.0 / j0_n);


   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   //real_t jp_ex = mp * pow(( (cp - cpmax) / (cp * ce)),0.5);
   //real_t jn_ex = mn * pow(( (cn - cpmax) / (cn * ce)),0.5);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t voltage = Up - Un  + eta_p - eta_n;

   // Temporary printing.
   if (Mpi::Root())
   {
      std::cout << "Up = " << Up << std::endl;
      std::cout << "Un = " << Un << std::endl;

      std::cout << "T = " << T << std::endl;
      std::cout << "I = " << I << std::endl;
      std::cout << "AP = " << AP << std::endl;
      std::cout << "AN = " << AN << std::endl;
      std::cout << "LPE = " << LPE << std::endl;
      std::cout << "LNE = " << LNE << std::endl;
      std::cout << "j0_p = " << j0_p << std::endl;
      std::cout << "j0_n = " << j0_n << std::endl;

      std::cout << "Voltage = " << voltage << std::endl;

      // Print file headings first time function is called.
      static bool writeFileHeadings = true;
      if (writeFileHeadings) {
         file << "t" << ", \t"
           << "voltage" << ", \t"
           << "cp" << ", \t"
           << "cn" << ", \t"
           << "eta_p" << ", \t"
           << "eta_n" << ", \t"
           << "Up" << ", \t"
           << "Un"
           << std::endl;

         writeFileHeadings = false;
      }

      // Print data to file.
      file << t << ", \t"
           << voltage << ", \t"
           << cp << ", \t" 
           << cn << ", \t" 
           << eta_p << ", \t"
           << eta_n << ", \t"
           << Up << ", \t"
           << Un
           << std::endl;
   }
}

void P2DOperator::GetParticleDofs(Array<int> & particle_dofs, Array<int> & particle_offsets)
{
   std::set<int> sep_gdofs_set;
   for (int e = 0; e < x_fespace->GetNE(); e++)
      if (x_fespace->GetAttribute(e) == SEP)
      {
         Array<int> dofs;
         x_fespace->GetElementDofs(e, dofs);
         for (int d: dofs)
            sep_gdofs_set.insert(x_fespace->GetGlobalTDofNumber(d));
      }

   unsigned max_sep_gdofs = NSEP * x_fespace->GetElementOrder(0) + 1;
   Array<int> sep_gdofs(max_sep_gdofs); sep_gdofs = -1;
   std::copy(sep_gdofs_set.begin(), sep_gdofs_set.end(), sep_gdofs.begin());

   Array<int> all_sep_gdofs(max_sep_gdofs * Mpi::WorldSize());
   MPI_Allgather(sep_gdofs.GetData(), max_sep_gdofs, MPI_INT,
                 all_sep_gdofs.GetData(), max_sep_gdofs, MPI_INT, MPI_COMM_WORLD);

   Array<int> boundary_dofs;
   x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (int dof = 0; dof < x_fespace->GetNDofs(); dof++)
   {
      int ltdof = x_fespace->GetLocalTDofNumber(dof);
      int gtdof = x_fespace->GetGlobalTDofNumber(dof);
      if (ltdof != -1 && boundary_dofs.Find(ltdof) == -1 && all_sep_gdofs.Find(gtdof) == -1)
         particle_dofs.Append(dof);
   }

   int my_particles = particle_dofs.Size();
   particle_offsets.SetSize(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

   particle_offsets.PartialSum();
   assert(particle_offsets[Mpi::WorldSize() - 1] == NPAR);
}
