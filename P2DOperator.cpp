#include "P2DOperator.hpp"
#include "utils.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace,
                         const unsigned &ndofs, BlockVector &x)
   : TimeDependentOperator(ndofs, (real_t) 0.0), x_fespace(x_fespace), r_fespace(r_fespace),
     A(NULL), current_dt(0.0), Solver(x_fespace->GetComm()), file("data.txt")
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
      GetParticleLTDofs(particle_dofs, particle_offsets);
      for (unsigned p = 0; p < NPAR; p++)
      {
         unsigned rank = std::distance(particle_offsets.begin(),
            std::find_if(particle_offsets.begin(), particle_offsets.end(), [&](int i){ return p < i; }));
         unsigned offset = Mpi::WorldRank() > 0 ? particle_offsets[Mpi::WorldRank() - 1] : 0;
         bool mine = p >= offset && p < offset + particle_dofs.Size();
         int ltdof = mine ? particle_dofs[p - offset] : -1;
         sc.Append(new SolidConcentration(*r_fespace[p], p, rank, ltdof));
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

   ep->Update(x);
   ec->Update(x);
   sp->Update(x);
   for (unsigned p = 0; p < NPAR; p++)
      sc[p]->Update(x);
}

real_t P2DOperator::ComputeExternalCurrent(const BlockVector &x)
{
   ParGridFunction cs_gf(x_fespace);
   cs_gf = 0;

   ParGridFunction ec_gf(x_fespace);
   ec_gf.SetFromTrueDofs(x.GetBlock(EC));

   for (unsigned p = 0; p < NPAR; p++)
   {
      real_t csurf = sc[p]->SurfaceConcentration(x);
      if (!isnan(csurf))
         cs_gf(sc[p]->GetParticleLTDof()) = csurf;
   }

   ExternalCurrentCoefficient spme_coeff(cs_gf, ec_gf);
   ParLinearForm sum(x_fespace);
   sum.AddDomainIntegrator(new DomainLFIntegrator(spme_coeff));
   sum.Assemble();

   real_t reduction_result = sum.Sum();
   MPI_Allreduce(MPI_IN_PLACE, &reduction_result, 1, MFEM_MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);

   return reduction_result / (LPE + LSEP + LNE);
}

void P2DOperator::ComputeVoltage(const BlockVector &x, real_t t, real_t dt)
{
   if (t == dt)
      file << "t" << ", "
           << "\t" << "voltage" << ", "
           << "\t" << "theta_p" << ", "
           << "\t" << "theta_n" << ", "
           << "\t" << "Up" << ", "
           << "\t" << "Un"
           << std::endl;

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

   //real_t R = 1.;//8.314;
   //real_t F = 1.;//96485;
   real_t T = 1.;//300.;
   //real_t Kp = 1.;
   //real_t Kn = 1.;
   real_t ce = 1.;//000.;
   //real_t ce0 = 1.;//000.;
   real_t cpmax = 1.;//30555;
   real_t cnmax = 1.;//51554;
   real_t Lp = 1./3.;
   real_t Ln = 1./3.;
   real_t Ap = 1.;
   real_t An = 1.;
   real_t I = 1.;
   real_t mp = 1.;
   real_t mn = 1.;

   real_t cp = csurf[0];   // Particle surface concentration at the positive electrode.
   real_t cn = csurf[1];   // Particle surface concentration at the negative electrode.

   //real_t jp0 = F * Kp * pow((ce/ce0) * (cp/cpmax) * (1 - (cp/cpmax)),0.5);
   //real_t jn0 = F * Kn * pow((ce/ce0) * (cn/cnmax) * (1 - (cn/cnmax)),0.5);

   // Definition from LIONSIMBA: https://doi.org/10.1149/2.0291607jes
   real_t theta_p = cp / cpmax;
   real_t theta_n = cn / cnmax;

   // Open Circuit Potential (no temperature dependence).
   // Definition from LIONSIMBA: https://doi.org/10.1149/2.0291607jes
   real_t Up_num = -4.656 + 88.669 * pow(theta_p,2) - 401.119 * pow(theta_p,4) +
                        342.909 * pow(theta_p,6) - 462.471 * pow(theta_p,8) + 433.434 * pow(theta_p,10);
   real_t Up_den = -1 + 18.933 * pow(theta_p,2) - 79.532 * pow(theta_p,4) +
                        37.311 * pow(theta_p,6) - 73.083 * pow(theta_p,8) + 95.96 * pow(theta_p,10);
   real_t Up = Up_num / Up_den;

   real_t Un = 0.7222 + 0.1387 * theta_n + 0.029 * pow(theta_n,0.5) - 0.0172 / theta_n +
                     0.0019 * pow(theta_n,-1.5) + 0.2808 * exp(0.9 - 15*theta_n) - 0.7984 * exp(0.4465 * theta_n - 0.4108);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t jp_ex = mp * pow(( (cp - cpmax) / (cp * ce)),0.5);
   real_t jn_ex = mn * pow(( (cn - cpmax) / (cn * ce)),0.5);

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t voltage = Up - Un  + 2 * T * (
                     asinh( I / (2 * Ap * Lp * jp_ex )) -
                     asinh( -I / (2 * An * Ln * jn_ex ) ) );

   // Temporary printing.
   if (Mpi::Root())
   {
      std::cout << "Up = " << Up << std::endl;
      std::cout << "Un = " << Un << std::endl;

      std::cout << "T = " << T << std::endl;
      std::cout << "I = " << I << std::endl;
      std::cout << "Ap = " << Ap << std::endl;
      std::cout << "An = " << An << std::endl;
      std::cout << "Lp = " << Lp << std::endl;
      std::cout << "Ln = " << Ln << std::endl;
      std::cout << "jp_ex = " << jp_ex << std::endl;
      std::cout << "jn_ex = " << jn_ex << std::endl;

      std::cout << "Voltage = " << voltage << std::endl;

      // Print data to file.
      file << t << ", "
           << "\t" << voltage << ", "
           << "\t" << theta_p << ", "
           << "\t" << theta_n << ", "
           << "\t" << Up << ", "
           << "\t" << Un
           << std::endl;
   }
}

void P2DOperator::GetParticleLTDofs(Array<int> & particle_dofs, Array<int> & particle_offsets)
{
   std::set<int> sep_global_dofs_set;
   for (int e = 0; e < x_fespace->GetNE(); e++)
      if (x_fespace->GetAttribute(e) == SEP)
      {
         Array<int> dofs;
         x_fespace->GetElementDofs(e, dofs);
         for (int d: dofs)
            sep_global_dofs_set.insert(x_fespace->GetGlobalTDofNumber(d));
      }

   unsigned max_sep_global_dofs = NSEP * x_fespace->GetElementOrder(0) + 1;
   Array<int> sep_global_dofs(max_sep_global_dofs); sep_global_dofs = -1;
   std::copy(sep_global_dofs_set.begin(), sep_global_dofs_set.end(), sep_global_dofs.begin());

   int all_sep_global_dofs[max_sep_global_dofs * Mpi::WorldSize()];
   MPI_Allgather(sep_global_dofs.GetData(), max_sep_global_dofs, MPI_INT,
                 all_sep_global_dofs, max_sep_global_dofs, MPI_INT, MPI_COMM_WORLD);

   sep_global_dofs_set = std::set<int>(all_sep_global_dofs, all_sep_global_dofs + max_sep_global_dofs * Mpi::WorldSize());

   std::set<int> particle_dofs_set;

   for (int d = 0; d < x_fespace->GetNDofs(); d++)
      if (sep_global_dofs_set.find(x_fespace->GetGlobalTDofNumber(d)) == sep_global_dofs_set.end())
         particle_dofs_set.insert(x_fespace->GetLocalTDofNumber(d));

   Array<int> boundary_dofs;
   x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (int d: boundary_dofs)
      particle_dofs_set.erase(d);
   particle_dofs_set.erase(-1);

   particle_dofs.SetSize(particle_dofs_set.size());
   std::copy(particle_dofs_set.begin(), particle_dofs_set.end(), particle_dofs.begin());

   int my_particles = particle_dofs.Size();
   particle_offsets.SetSize(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

   particle_offsets.PartialSum();
   assert(particle_offsets[Mpi::WorldSize() - 1] == NPAR);
}
