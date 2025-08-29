#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &x)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     Ac(NULL), Ap(NULL), current_dt(0.0), Solver(x_fespace->GetComm()), file("data.csv")
{
   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   //Solver.SetPreconditioner(Prec);

   block_trueOffsets.SetSize(NEQS + 1);
   potential_trueOffsets.SetSize(NMACROP + 1);
   concentration_trueOffsets.SetSize(NMACROC + NPAR + 1);

   block_trueOffsets[0] = 0;
   block_trueOffsets[EP + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[SP + 1] = x_fespace->GetTrueVSize();
   block_trueOffsets[EC + 1] = x_fespace->GetTrueVSize();

   potential_trueOffsets[0] = 0;
   potential_trueOffsets[EPP + 1] = x_fespace->GetTrueVSize();
   potential_trueOffsets[SPP + 1] = x_fespace->GetTrueVSize();

   concentration_trueOffsets[0] = 0;
   concentration_trueOffsets[ECC + 1] = x_fespace->GetTrueVSize();

   for (unsigned p = 0; p < NPAR; p++)
   {
      block_trueOffsets[SC + 1 + p] = r_fespace[p]->GetTrueVSize();
      concentration_trueOffsets[SCC + 1 + p] = r_fespace[p]->GetTrueVSize();
   }

   block_trueOffsets.PartialSum();
   potential_trueOffsets.PartialSum();
   concentration_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << NEQS << std::endl;
      std::cout << "Unknowns (rank 0): " << block_trueOffsets[NEQS] << std::endl;
   }

   x.Update(block_trueOffsets); x = 0.;
   bp.Update(potential_trueOffsets);
   bc.Update(concentration_trueOffsets);

   ep = new ElectrolytePotential(*x_fespace);
   ec = new ElectrolyteConcentration(*x_fespace);
   sp = new SolidPotential(*x_fespace);

   if (M == SPM || M == SPMe)
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

   Array<int> offsets({block_trueOffsets[P], block_trueOffsets[C], block_trueOffsets[NEQS]});
   BlockVector dx_dt_blocked(dx_dt, offsets);

   // assemble Ac and bc, create dxc_dt ref
   Ac->SetDiagonalBlock(ECC, Add(1, ec->GetM(), dt, ec->GetK()));
   bc.GetBlock(ECC) = ec->GetZ();
   for (unsigned p = 0; p < NPAR; p++)
   {
      Ac->SetDiagonalBlock(SCC + p, Add(1, sc[p]->GetM(), dt, sc[p]->GetK()));
      bc.GetBlock(SCC + p) = sc[p]->GetZ();
   }
   Vector & dxc_dt(dx_dt_blocked.GetBlock(1));

   // solve for dxc_dt (concentrations rate)
   Solver.SetOperator(*Ac);
   Solver.Mult(bc, dxc_dt);

   // assemble Ap and bp, create dxp_dt ref
   Ap->SetDiagonalBlock(EPP, new HypreParMatrix(ep->GetK()));
   Ap->SetDiagonalBlock(SPP, new HypreParMatrix(sp->GetK()));
   bp.GetBlock(EPP) = ep->GetZ();
   bp.GetBlock(SPP) = sp->GetZ();
   Vector & dxp_dt(dx_dt_blocked.GetBlock(0));

   // solve for dxp_dt (potentials rate)
   Solver.SetOperator(*Ap);
   Solver.Mult(bp, dxp_dt);
}

void P2DOperator::Update(const BlockVector &x, const real_t &dt)
{
   // rebuild Ap, Ac
   delete Ap;
   delete Ac;
   Ap = new BlockOperator(potential_trueOffsets);
   Ac = new BlockOperator(concentration_trueOffsets);
   Ap->owns_blocks = 1;
   Ac->owns_blocks = 1;

   ConstantCoefficient jx = ComputeReactionCurrent(x);
   ep->Update(x, jx, dt);
   sp->Update(x, jx, dt);

   if (M == SPM || M == SPMe)
   {
      ec->Update(x, ComputeReactionCurrent());
      for (unsigned p = 0; p < NPAR; p++)
         sc[p]->Update(x, ComputeReactionCurrent(sc[p]->GetParticleRegion()));
   }
   else if (M == P2D)
   {
      ParGridFunction j(x_fespace);
      j.ProjectCoefficient(jx);

      for (unsigned p = 0; p < NPAR; p++)
      {
         MPI_Request request;
         real_t jr = sc[p]->IsParticleOwned() ? j(sc[p]->GetParticleDof()) : 0;
         if (sc[p]->IsParticleOwned())
            MPI_Isend(&jr, 1, MFEM_MPI_REAL_T, sc[p]->GetSurfaceRank(), 1, MPI_COMM_WORLD, &request);

         if (sc[p]->IsSurfaceOwned())
            MPI_Recv(&jr, 1, MFEM_MPI_REAL_T, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         if (sc[p]->IsParticleOwned())
            MPI_Wait(&request, MPI_STATUS_IGNORE);

         sc[p]->Update(x, ConstantCoefficient(jr));
      }
   }
}

ConstantCoefficient P2DOperator::ComputeReactionCurrent(const Region &r)
{
   return ConstantCoefficient(ComputeReactionCurrent()(r));
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

real_t P2DOperator::GetSurfaceConcentration(const Region &r, const BlockVector &x)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant surface concentration, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == PE || r == NE, "Cannot get constant surface concentration, only positive (PE) and negative electrodes (NE) are supported.");

   real_t csurf = r == PE ? CP0 : CN0;
   for (unsigned p = 0; p < NPAR; p++)
      if (sc[p]->GetParticleRegion() == r)
         csurf += sc[p]->SurfaceConcentration(x);

   return csurf;
}

ParGridFunction P2DOperator::GetSurfaceConcentration(const BlockVector &x)
{
   ParGridFunction sc_gf(x_fespace);
   sc_gf = 0;

   for (unsigned p = 0; p < NPAR; p++)
      if (sc[p]->IsParticleOwned())
      {
         real_t csurf0 = sc[p]->GetParticleRegion() == PE ? CP0 : CN0;
         sc_gf(sc[p]->GetParticleDof()) = csurf0 + sc[p]->SurfaceConcentration(x);
      }

   // Apply prolongation after restriction. Might be unnecessary, but guarantees
   // all processors have the right information for all their local dofs.
   sc_gf.SetFromTrueVector();

   return sc_gf;
}

real_t P2DOperator::ComputeExchangeCurrent(const Region &r, const BlockVector &x)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant exchange current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == PE || r == NE, "Cannot calculate constant exchange current, only positive (PE) and negative electrodes (NE) are supported.");

   real_t sc = GetSurfaceConcentration(r, x);
   real_t k = r == PE ? KP : KN;

   if (M == SPM)
      return ExchangeCurrentCoefficient(k, sc, CE0).Eval();
   else
   {
      ParGridFunction ec_gf(x_fespace);
      ec_gf.SetFromTrueDofs(x.GetBlock(EC));

      ExchangeCurrentCoefficient jex(k, sc, ec_gf);

      Array<int> markers(x_fespace->GetParMesh()->attributes.Max());
      markers = 0; markers[r-1] = 1;
      ParLinearForm sum(x_fespace);
      sum.AddDomainIntegrator(new DomainLFIntegrator(jex), markers);
      sum.Assemble();

      real_t l = r == PE ? LPE : LNE;
      return x_fespace->GetParMesh()->ReduceInt(sum.Sum()) / l;
   }
}

ExchangeCurrentCoefficient P2DOperator::ComputeExchangeCurrent(const BlockVector &x)
{
   ParGridFunction sc_gf = GetSurfaceConcentration(x);
   ParGridFunction ec_gf(x_fespace);
   ec_gf.SetFromTrueDofs(x.GetBlock(EC));

   return ExchangeCurrentCoefficient(1, sc_gf, ec_gf);
}

real_t P2DOperator::ComputeOpenCircuitPotential(const Region &r, const real_t &x)
{
   MFEM_ASSERT(r == PE || r == NE, "Cannot calculate open circuit potential, only positive (PE) and negative electrodes (NE) are supported.")
   return r == PE ? UP(x) : UN(x);
}

void P2DOperator::ComputeVoltage(const BlockVector &x, real_t t, real_t dt)
{
   real_t Up = ComputeOpenCircuitPotential(PE, GetSurfaceConcentration(PE, x));
   real_t Un = ComputeOpenCircuitPotential(NE, GetSurfaceConcentration(NE, x));

   real_t jp = ComputeReactionCurrent(PE).constant;
   real_t jn = ComputeReactionCurrent(NE).constant;

   real_t j0_p =  ComputeExchangeCurrent(PE, x);
   real_t j0_n =  ComputeExchangeCurrent(NE, x);

   real_t eta_p = 2 * T * asinh(jp / 2.0 / j0_p);
   real_t eta_n = 2 * T * asinh(jn / 2.0 / j0_n);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t voltage = Up - Un + (eta_p - eta_n) * phi_scale;

   // Temporary printing (a bit miraculous this is working even with GetSurfaceConcentration but causes MPI error on termination)
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
           << GetSurfaceConcentration(PE, x) << ", \t"
           << GetSurfaceConcentration(NE, x) << ", \t"
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
   MFEM_ASSERT(particle_offsets[Mpi::WorldSize() - 1] == NPAR, "Failed to distribute particles across processors.");
}
