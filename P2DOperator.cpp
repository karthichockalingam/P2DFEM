#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &x, const real_t & dt)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     Ac(NULL), Ap(NULL), _x(x), _dt(dt), Solver(x_fespace->GetComm()), file("data.csv")
{
   const real_t rel_tol = 1e-8;

   Solver.iterative_mode = false;
   Solver.SetRelTol(rel_tol);
   Solver.SetAbsTol(0.0);
   Solver.SetMaxIter(100);
   Solver.SetPrintLevel(0);
   //Solver.SetPreconditioner(Prec);

   block_offsets.SetSize(NMACRO + 1 + 1);
   block_trueOffsets.SetSize(NEQS + 1);
   potential_trueOffsets.SetSize(NMACROP + 1);
   concentration_trueOffsets.SetSize(NMACROC + NPAR + 1);

   block_offsets[0] = 0;
   block_offsets[EP + 1] = x_fespace->GetVSize();
   block_offsets[SP + 1] = x_fespace->GetVSize();
   block_offsets[EC + 1] = x_fespace->GetVSize();
   block_offsets[SC + 1] = x_fespace->GetVSize();

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

   block_offsets.PartialSum();
   block_trueOffsets.PartialSum();
   potential_trueOffsets.PartialSum();
   concentration_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << NEQS << std::endl;
      std::cout << "Unknowns (rank 0): " << block_trueOffsets[NEQS] << std::endl;
   }

   // Set offsets for full dof vector
   _l.Update(block_offsets); _l = 0.;

   // Initialise gridfunctions to use the appropriate section of the full dof vector _l
   _ep_gf = new ParGridFunction(x_fespace, _l, block_offsets[EP]);
   _ec_gf = new ParGridFunction(x_fespace, _l, block_offsets[SP]);
   _sp_gf = new ParGridFunction(x_fespace, _l, block_offsets[EC]);
   _sc_gf = new ParGridFunction(x_fespace, _l, block_offsets[SC]);

   // Set offsets for solution and rhs (potential and concentration) true vectors
   _x.Update(block_trueOffsets); _x = 0.;
   bp.Update(potential_trueOffsets);
   bc.Update(concentration_trueOffsets);

   // Construct equation ojects, first the 3 macro equations, then the NPAR micro eqs
   ep = new ElectrolytePotential(*x_fespace);
   sp = new SolidPotential(*x_fespace);
   ec = new ElectrolyteConcentration(*x_fespace);

   if (M == SPM || M == SPMe)
   {
      sc.Append(new SolidConcentration(*r_fespace[0], 0, 0, -1, NE));
      sc.Append(new SolidConcentration(*r_fespace[1], 1, 0, -1, PE));
   }
   else
   {
      Array<int> particle_dofs, particle_offsets;
      Array<Region> particle_regions;
      GetParticleDofs(particle_dofs, particle_regions, particle_offsets);

      for (unsigned p = 0; p < NPAR; p++)
      {
         auto rank_iter = std::upper_bound(particle_offsets.begin(), particle_offsets.end(), p);
         unsigned rank = std::distance(particle_offsets.begin(), rank_iter) - 1;
         bool owned = rank == Mpi::WorldRank();

         unsigned offset = particle_offsets[Mpi::WorldRank()];
         int dof = owned ? particle_dofs[p - offset] : -1;
         Region region = owned ? particle_regions[p - offset] : UNKNOWN;

         sc.Append(new SolidConcentration(*r_fespace[p], p, rank, dof, region));
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

   if (M == P2D)
   {
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
}

void P2DOperator::SetGridFunctionsFromTrueVectors()
{
   _ep_gf->SetFromTrueDofs(_x.GetBlock(EP));
   _sp_gf->SetFromTrueDofs(_x.GetBlock(SP));
   _ec_gf->SetFromTrueDofs(_x.GetBlock(EC));
   SetSurfaceConcentration();
}

void P2DOperator::Update()
{
   // rebuild Ap, Ac
   delete Ap;
   delete Ac;
   Ap = new BlockOperator(potential_trueOffsets);
   Ac = new BlockOperator(concentration_trueOffsets);
   Ap->owns_blocks = 1;
   Ac->owns_blocks = 1;

   if (M == SPM || M == SPMe)
   {
      ec->Update(_x, ComputeReactionCurrent());
      for (unsigned p = 0; p < NPAR; p++)
         sc[p]->Update(_x, ConstantCoefficient(ComputeReactionCurrent(sc[p]->GetParticleRegion())));
   }
   else if (M == P2D)
   {
      ReactionCurrentCoefficient jx = ComputeReactionCurrent();
      ep->Update(_x, jx, _dt);
      sp->Update(_x, jx, _dt);
      ec->Update(_x, jx, _dt);
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

         sc[p]->Update(_x, ConstantCoefficient(jr));
      }
   }
}

//
// Surface Concentration
//

real_t P2DOperator::GetSurfaceConcentration(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant surface concentration, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant surface concentration, only negative (NE) and positive electrodes (PE) are supported.");

   real_t csurf = r == NE ? CN0 : CP0;
   for (unsigned p = 0; p < NPAR; p++)
      if (sc[p]->GetParticleRegion() == r)
         csurf += sc[p]->SurfaceConcentration(_x);

   return csurf;
}

void P2DOperator::SetSurfaceConcentration()
{
   if (M == SPM || M == SPMe)
      return;

   for (unsigned p = 0; p < NPAR; p++)
      if (sc[p]->IsParticleOwned())
      {
         real_t csurf0 = sc[p]->GetParticleRegion() == NE ? CN0 : CP0;
         (*_sc_gf)(sc[p]->GetParticleDof()) = csurf0 + sc[p]->SurfaceConcentration(_x);
      }

   // Apply prolongation after restriction. Might be unnecessary, but guarantees
   // all processors have the right information for all their local dofs.
   _sc_gf->SetFromTrueVector();
}

//
// Reaction Current
//

real_t P2DOperator::ComputeReactionCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant reaction current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant reaction current, only negative (NE) and positive electrodes (PE) are supported.")
   return ComputeReactionCurrent().Eval()(r);
}

ReactionCurrentCoefficient P2DOperator::ComputeReactionCurrent()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         return ReactionCurrentCoefficient();
      case P2D:
         return ReactionCurrentCoefficient(T, ComputeExchangeCurrent(), ComputeOverPotential());
   }
   MFEM_ASSERT(false, "Unreachable.");
}

//
// Exchange Current
//

real_t P2DOperator::ComputeExchangeCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant exchange current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant exchange current, only negative (NE) and positive electrodes (PE) are supported.");
   return ComputeExchangeCurrent().Eval()(r);
}

ExchangeCurrentCoefficient P2DOperator::ComputeExchangeCurrent()
{
   switch (M)
   {
      case SPM:
         return ExchangeCurrentCoefficient(KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), CE0);
      case SPMe:
         return ExchangeCurrentCoefficient(KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), *_ec_gf);
      case P2D:
         return ExchangeCurrentCoefficient(KN, KP, *_sc_gf, *_ec_gf);
   }
   MFEM_ASSERT(false, "Unreachable.");
}

//
// Open Circuit Potential
//

real_t P2DOperator::ComputeOpenCircuitPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant open circuit potential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant open circuit potential, only negative (NE) and positive electrodes (PE) are supported.")
   return ComputeOpenCircuitPotential().Eval()(r);
}

OpenCircuitPotentialCoefficient P2DOperator::ComputeOpenCircuitPotential()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         return OpenCircuitPotentialCoefficient(UN, UP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE));
      case P2D:
         return OpenCircuitPotentialCoefficient(UN, UP, *_sc_gf);
   }
   MFEM_ASSERT(false, "Unreachable.");
}

//
// OverPotential
//

real_t P2DOperator::ComputeOverPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant  overpotential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant overpotential, only negative (NE) and positive electrodes (PE) are supported.")
   return ComputeOverPotential().Eval()(r);
}

OverPotentialCoefficient P2DOperator::ComputeOverPotential()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         return OverPotentialCoefficient(T, ComputeExchangeCurrent());
      case P2D:
         return OverPotentialCoefficient(*_sp_gf, *_ep_gf, ComputeOpenCircuitPotential());
   }
   MFEM_ASSERT(false, "Unreachable.");
}

//
// Voltage
//

void P2DOperator::ComputeVoltage(const BlockVector &x, real_t t, real_t dt)
{
   real_t Un = ComputeOpenCircuitPotential(NE);
   real_t Up = ComputeOpenCircuitPotential(PE);

   real_t eta_n = ComputeOverPotential(NE);
   real_t eta_p = ComputeOverPotential(PE);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t voltage = Up - Un + (eta_p - eta_n) * phi_scale;

   real_t scn = GetSurfaceConcentration(NE);
   real_t scp = GetSurfaceConcentration(PE);
   if (Mpi::Root())
   {
      std::cout << "[Rank " << Mpi::WorldRank() << "]"
                << " Voltage = " << voltage << std::endl;

      // Print file headings first time function is called.
      static bool writeFileHeadings = true;
      if (writeFileHeadings) {
         file << "t" << ", \t"
           << "cn" << ", \t"
           << "cp" << ", \t"
           << "voltage"
           << std::endl;

         writeFileHeadings = false;
      }

      // Print data to file.
      file << t << ", \t"
           << scn << ", \t"
           << scp << ", \t"
           << voltage
           << std::endl;
   }

   sc[1]->DebuggingCheck(_x);
}

void P2DOperator::GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets)
{
   std::set<std::pair<int, Region>> electrode_dofs_set;
   std::set<int> sep_gdofs_set;
   for (int e = 0; e < x_fespace->GetNE(); e++)
   {
      Array<int> dofs;
      x_fespace->GetElementDofs(e, dofs);
      switch (Region r = Region(x_fespace->GetAttribute(e)))
      {
         case NE:
         case PE:
            for (int d: dofs)
               electrode_dofs_set.insert({d, r});
            break;
         case SEP:
            for (int d: dofs)
               sep_gdofs_set.insert(x_fespace->GetGlobalTDofNumber(d));
            break;
      }
   }

   unsigned max_sep_gdofs = NSEP * x_fespace->FEColl()->GetOrder() + 1;
   Array<int> sep_gdofs(max_sep_gdofs); sep_gdofs = -1;
   std::copy(sep_gdofs_set.begin(), sep_gdofs_set.end(), sep_gdofs.begin());

   Array<int> all_sep_gdofs(max_sep_gdofs * Mpi::WorldSize());
   MPI_Allgather(sep_gdofs.GetData(), max_sep_gdofs, MPI_INT,
                 all_sep_gdofs.GetData(), max_sep_gdofs, MPI_INT, MPI_COMM_WORLD);

   Array<int> boundary_dofs;
   x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (auto [dof, region]: electrode_dofs_set)
   {
      int ltdof = x_fespace->GetLocalTDofNumber(dof);
      int gtdof = x_fespace->GetGlobalTDofNumber(dof);
      if (ltdof != -1 && boundary_dofs.Find(ltdof) == -1 && all_sep_gdofs.Find(gtdof) == -1)
      {
         particle_dofs.Append(dof);
         particle_regions.Append(region);
      }
   }

   int my_particles = particle_dofs.Size();
   particle_offsets.SetSize(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

   particle_offsets.Prepend(0);
   particle_offsets.PartialSum();
   MFEM_ASSERT(particle_offsets[Mpi::WorldSize()] == NPAR, "Failed to distribute particles across processors.");
}
