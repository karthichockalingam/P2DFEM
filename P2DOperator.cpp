#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, real_t dt, BlockVector &x)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     A(NULL), current_dt(dt), Solver(x_fespace->GetComm()), file("data.csv")
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

   ConstantCoefficient jx = ComputeReactionCurrent(x);
   ep->Update(x, jx, current_dt);
   sp->Update(x, jx, current_dt);

   if (M == SPM || M == SPMe)
   {
      ec->Update(x, ComputeReactionCurrent(), current_dt);
      for (unsigned p = 0; p < NPAR; p++)
         sc[p]->Update(x, ComputeReactionCurrent(sc[p]->GetParticleRegion()), current_dt);
   }
   else if (M == P2D)
   {
      ParGridFunction j(x_fespace);
      j.ProjectCoefficient(jx);

      for (unsigned p = 0; p < NPAR; p++)
      {
         // incomplete, needs comm
         ConstantCoefficient jr(j(sc[p]->GetParticleDof()));
         sc[p]->Update(x, jr, current_dt);
      }
   }
}

ConstantCoefficient P2DOperator::ComputeReactionCurrent(const Region &r)
{
   return ConstantCoefficient(ComputeReactionCurrent()(r + 1));
}

PWConstCoefficient P2DOperator::ComputeReactionCurrent()
{
   Vector c({/* PE */ - I / AP / LPE, /* SEP */ 0., /* NE */ + I / AN / LNE});
   return PWConstCoefficient(c);
}

ConstantCoefficient P2DOperator::ComputeReactionCurrent(const BlockVector &x)
{
   return ConstantCoefficient();
}

real_t P2DOperator::ComputeSurfaceConcentration(const Region &r, const BlockVector &x)
{
   real_t sc0 = r == PE ? CP0 : r == NE ? CN0 : 0;
   for (unsigned p = 0; p < NPAR; p++)
      if (sc[p]->GetParticleRegion() == r)
         return sc0 + sc[p]->SurfaceConcentration(x);

   mfem_error("Cannot calculate constant surface concentration for the given region. Only positive (PE) and negative electrodes (NE) are supported.");
}

ParGridFunction P2DOperator::ComputeSurfaceConcentration(const BlockVector &x)
{
   ParGridFunction sc_gf(x_fespace);
   sc_gf = 0;

   for (unsigned p = 0; p < NPAR; p++)
   {
      Region r = sc[p]->GetParticleRegion();
      real_t sc0 = r == PE ? CP0 : r == NE ? CN0 : 0;
      real_t csurf = sc[p]->SurfaceConcentration(x);
      if (sc[p]->IsParticleOwned())
         sc_gf(sc[p]->GetParticleDof()) = sc0 + csurf;
   }
   // Apply prolongation after restriction. Might be unnecessary, but guarantees
   // all processors have the right information for all their local dofs.
   sc_gf.SetFromTrueVector();

   return sc_gf;
}

real_t P2DOperator::ComputeExchangeCurrent(const Region &r, const BlockVector &x)
{
   if (M == SPM)
   {
      real_t sc = ComputeSurfaceConcentration(r, x);
      real_t k = r == PE ? KP : r == NE ? KN : 0;
      return k * sqrt(sc * CE0 * abs(1.0 - sc));
   }
   else if (M == SPMe)
   {
      ExchangeCurrentCoefficient jex = ComputeExchangeCurrent(x);
      Array<int> markers(x_fespace->GetParMesh()->attributes.Max());
      markers = 0; markers[r] = 1;
      ParLinearForm sum(x_fespace);
      sum.AddDomainIntegrator(new DomainLFIntegrator(jex), markers);
      sum.Assemble();

      real_t reduction_result = sum.Sum();
      MPI_Allreduce(MPI_IN_PLACE, &reduction_result, 1, MFEM_MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);

      real_t l = r == PE ? LPE : r == NE ? LNE : 0;
      return reduction_result / l;
   }
   else
      mfem_error("Cannot calculate constant exchange current for the given method. Only SPM and SPMe are supported.");
}

ExchangeCurrentCoefficient P2DOperator::ComputeExchangeCurrent(const BlockVector &x)
{
   ParGridFunction sc_gf = ComputeSurfaceConcentration(x);
   ParGridFunction ec_gf(x_fespace);
   ec_gf.SetFromTrueDofs(x.GetBlock(EC));

   return ExchangeCurrentCoefficient(sc_gf, ec_gf);
}

real_t P2DOperator::ComputeOpenCircuitPotential(const Region &r, const real_t &x)
{
   if (r == PE)
      return -0.8090*x + 4.4875 - 0.0428*tanh(18.5138*(x - 0.5542))
             - 17.7326*tanh(15.7890*(x - 0.3117)) + 17.5842*tanh(15.9308*(x - 0.3120));
   else if (r == NE)
      return 1.97938*exp(-39.3631*x) + 0.2482 - 0.0909*tanh(29.8538*(x - 0.1234))
             - 0.04478*tanh(14.9159*(x - 0.2769)) - 0.0205*tanh(30.4444*(x - 0.6103));
   else
      mfem_error("Cannot calculate open circuit potential for the given region. Only positive (PE) and negative electrodes (NE) are supported.");
}

void P2DOperator::ComputeVoltage(const BlockVector &x, real_t t, real_t dt)
{
   real_t Up = ComputeOpenCircuitPotential(PE, ComputeSurfaceConcentration(PE, x));
   real_t Un = ComputeOpenCircuitPotential(NE, ComputeSurfaceConcentration(NE, x));

   real_t jp = ComputeReactionCurrent(PE).constant;
   real_t jn = ComputeReactionCurrent(NE).constant;

   real_t j0_p =  ComputeExchangeCurrent(PE, x);
   real_t j0_n =  ComputeExchangeCurrent(NE, x);

   real_t eta_p = 2 * T * asinh(jp / 2.0 / j0_p);
   real_t eta_n = 2 * T * asinh(jn / 2.0 / j0_n);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t voltage = Up - Un + (eta_p - eta_n) * phi_scale;

   // Temporary printing.
   if (Mpi::Root())
   {
      std::cout << "[Rank " << Mpi::WorldRank() << "]"
                  << " Voltage = " << voltage << std::endl;

      // Print file headings first time function is called.
      static bool writeFileHeadings = true;
      if (writeFileHeadings) {
         file << "t" << ", \t"
           << "cp" << ", \t"
           << "cn" << ", \t"
           << "voltage"
           << std::endl;

         writeFileHeadings = false;
      }

      // Print data to file.
      file << t << ", \t"
           << ComputeSurfaceConcentration(PE, x) << ", \t"
           << ComputeSurfaceConcentration(NE, x) << ", \t"
           << voltage
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
