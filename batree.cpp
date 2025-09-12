// Batree - An MFEM-based SPM, SPMe and P2D solver by the Hartree Centre
//
// Compile with: make batree
//
// Sample runs:  mpirun -np 4 batree -m SPM
//               mpirun -np 4 batree -m SPMe
//
// Description:  Under active development. No explicit time integration methods
//               are supported at this time. Use -m or --method to select from
//               the three electrochemical models.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "P2DOperator.hpp"

using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   std::string model = "SPM";
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int order = 1;
   int ode_solver_type = 1;
   real_t t_final = 3600.0;
   real_t dt = 1.0;
   bool visualization = true;
   int vis_steps = 5;

   int precision = 8;
   std::cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&model, "-m", "--model",
                  "Electrochemical model: SPM, SPMe, or P2D.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                  "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable ParaView visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }

   if (Mpi::Root())
      args.PrintOptions(std::cout);

   // 4. Define the ODE solver used for time integration. Several implicit
   //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
   //    explicit Runge-Kutta methods are available.
   ODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver = new BackwardEulerSolver; break;
      case 2:  ode_solver = new SDIRK23Solver(2); break;
      case 3:  ode_solver = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver = new ForwardEulerSolver; break;
      case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver = new RK3SSPSolver; break;
      case 14: ode_solver = new RK4Solver; break;
      case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      default:
      std::cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
      return 1;
   }

   // Initialise grid and layout properties dependent on the electrochemical model and FE order
   init_params(model, order);

   // 3. Read the serial mesh from the given mesh file on all processors. We can
   //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
   //    with the same code.
   Mesh x_smesh = Mesh::MakeCartesian1D(NX);
   for (unsigned i = 0; i < NX; i++)
      x_smesh.SetAttribute(i, i < NNE ? NE : i < NNE + NSEP ? SEP : PE);
   x_smesh.SetAttributes();

   Mesh r_smesh[NPAR];
   for (unsigned p = 0; p < NPAR; p++)
      r_smesh[p] = Mesh::MakeCartesian1D(NR);

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *x_pmesh = new ParMesh(MPI_COMM_WORLD, x_smesh);
   x_smesh.Clear(); // the serial mesh is no longer needed
   ParMesh *r_pmesh[NPAR];
   for (unsigned p = 0; p < NPAR; p++)
   {
      r_pmesh[p] = new ParMesh(MPI_COMM_WORLD, r_smesh[p]);
      r_smesh[p].Clear(); // the serial mesh is no longer needed
   }

   // 7. Define the vector finite element space representing the current and the
   //    initial temperature, u_ref.
   H1_FECollection fe_coll(order, /*dim*/ 1);
   ParFiniteElementSpace * x_fespace = new ParFiniteElementSpace(x_pmesh, &fe_coll);
   Array<ParFiniteElementSpace *> r_fespace(NPAR);
   for (unsigned p = 0; p < NPAR; p++)
      r_fespace[p] = new ParFiniteElementSpace(r_pmesh[p], &fe_coll);

   // 8. Get the total number of dofs in the system (including boundaries)
   {
      HYPRE_BigInt fe_size_global = NMACRO * x_fespace->GlobalTrueVSize();
      for (unsigned p = 0; p < NPAR; p++)
         fe_size_global += r_fespace[p]->GlobalTrueVSize();

      if (Mpi::Root())
         std::cout << "Unknowns (total): " << fe_size_global << std::endl;
   }

   // 8.5 Get the total number of dofs _owned_ by this processor
   HYPRE_BigInt fe_size_owned = NMACRO * x_fespace->GetTrueVSize();
   for (unsigned p = 0; p < NPAR; p++)
      fe_size_owned += r_fespace[p]->GetTrueVSize();

   // 9. Initialize the conduction operator and the VisIt visualization.
   BlockVector x;
   P2DOperator oper(x_fespace, r_fespace, fe_size_owned, x, dt);

   // 10. Perform time-integration (looping over the time iterations, ti, with a
   //     time-step dt).
   ode_solver->Init(oper);
   real_t t = 0.0;

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      last_step = t + dt >= t_final - dt/2;

      ode_solver->Step(x, t, dt);
      oper.SetGridFunctionsFromTrueVectors();
      oper.ComputeVoltage(x, t, dt);

      if (last_step || (ti % vis_steps) == 0)
      {
         if (Mpi::Root())
            std::cout << "step " << ti << ", t = " << t << std::endl;

         // TODO: Stop sim at cutoff voltage
         if (last_step || visualization)
         {
         }
      }
   }

   // 11. Free the used memory.
   delete ode_solver;
   delete x_pmesh;
   for (unsigned p = 0; p < NPAR; p++)
      delete r_pmesh[p];

   return 0;
}
